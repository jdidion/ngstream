from abc import ABCMeta, abstractmethod
import base64
import io
import logging
import os
from pathlib import Path
from queue import Queue, Empty
import re
from threading import Thread
import time
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence, Union, cast
from urllib.parse import ParseResult, parse_qs, urlencode, urlunparse

import pysam
import requests
from xphyle import parse_url

from ngstream.api import Record


CONTENT_LENGTH = "Content-Length"
CONTENT_TYPE = "Content-Type"
HTSGET_KEY = "htsget"
URL_KEY = "url"
HEADER_CLASS = "header"

CONTENT_TYPE_RE = re.compile(r"application/vnd\.ga4gh\.htsget\.v(.*)+json.*")
REGION_RE = re.compile(r"([^:]*)(?::(.*?)-(.*))?")

DEFAULT_CHUNK_SIZE = 2 ** 16

BGZF_EOF_MARKER = bytearray(
    [
        0x1F,
        0x8B,
        0x08,
        0x04,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0xFF,
        0x06,
        0x00,
        0x42,
        0x43,
        0x02,
        0x00,
        0x1B,
        0x00,
        0x03,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
    ]
)
"""These 28 bytes signal the end of a BGZF block. We need to add them to the data
streamed from 'header' class URLs in order to create a valid BAM file that can be
parsed by pysam.
"""


class DownloaderStateError(Exception):
    """
    The Downloader is in an invalid state.
    """


class ContentLengthMismatchError(Exception):
    """
    The length of the downloaded content is not the same as the length reported in
    the header.
    """


class SamRecord(Record):
    @property
    @abstractmethod
    def is_paired(self) -> bool:
        pass

    @property
    @abstractmethod
    def is_read1(self) -> bool:
        pass

    @property
    @abstractmethod
    def is_read2(self) -> bool:
        pass


class HtsgetDownloader(Thread, metaclass=ABCMeta):
    """
    Thread that downloads BAM data from URLs, passes the data to 'samtools view' to
    convert to SAM (which runs in a second thread), and iterates over output records.
    """

    FINISH_SIGNAL = 0
    TERMINATE_SIGNAL = 1

    def __init__(self, timeout: int = 10):
        super().__init__()
        self.timeout = timeout
        self.daemon = True
        self._url_queue = Queue()

    # external interface

    def download_once(
        self, url: Union[str, ParseResult], *ticket_args, **ticket_kwargs
    ) -> None:
        ticket = get_ticket(url, *ticket_args, **ticket_kwargs)
        self.download_urls_once(ticket["urls"])

    def download_urls_once(self, url_objects: Sequence[dict]) -> None:
        """
        Starts the downloader, download a single list of URLs, and finish.

        Args:
            url_objects:
        """
        self.start()
        try:
            self.download_urls(url_objects)
            self.finish(now=False)
        except:
            self.finish(now=True)

    def download(
        self, url: Union[str, ParseResult], *ticket_args, **ticket_kwargs
    ) -> None:
        ticket = get_ticket(url, *ticket_args, **ticket_kwargs)
        self.download_urls(ticket["urls"])

    def download_urls(self, url_objects: Sequence[dict]) -> None:
        """
        Adds a URL to the download queue.
        """
        if not self.is_alive():
            raise DownloaderStateError("Downloader not started")

        for url_object in url_objects:
            self._url_queue.put(url_object)

        self._url_queue.put(HtsgetDownloader.FINISH_SIGNAL)

    def finish(self, now: bool = False) -> None:
        """
        Stops the Downloader.

        Args:
            now: Whether to force the downloader to stop. If True, the URL queue is
                first emptied. Then a TERMINATE_SIGNAL is added to the queue. This
                causes the thread to kill the SAM->BAM process and end immediately.
                Otherwise this blocks until the Downloader has finished with all the
                URLs currently in its queue.

        Raises:
             RuntimeError if `now` is True and the thread fails to terminate.
        """
        if now:
            # Empty out the queue
            while True:
                try:
                    self._url_queue.get_nowait()
                except Empty:
                    break

            self._url_queue.put(HtsgetDownloader.TERMINATE_SIGNAL)

            try:
                self.join(timeout=self.timeout)
            except TimeoutError:
                pass

            # Raise an error if we didn't die
            if self.is_alive():
                raise RuntimeError("Downloader did not terminate")
        else:
            self._url_queue.put(HtsgetDownloader.FINISH_SIGNAL)
            self.join()

    # thread interface

    def run(self):
        signal = HtsgetDownloader.FINISH_SIGNAL

        try:
            self._init_thread()

            while True:
                try:
                    url_object = self._url_queue.get_nowait()

                    if isinstance(url_object, int):
                        signal = url_object
                        break

                    url = parse_url(url_object[URL_KEY])

                    if url is None:
                        raise ValueError(f"Invalid URL: {url_object['url']}")

                    handle_htsget_url(url, url_object, self.timeout, self._write)
                except Empty:
                    time.sleep(1)
        finally:
            if signal is HtsgetDownloader.FINISH_SIGNAL:
                self._finish_thread()
            else:
                self._terminate_thread()

    @abstractmethod
    def _write(self, data: bytes) -> None:
        pass

    @abstractmethod
    def _init_thread(self):
        pass

    @abstractmethod
    def _finish_thread(self) -> None:
        pass

    def _terminate_thread(self) -> None:
        self._finish_thread()


class BaseSamHtsgetDownloader(HtsgetDownloader, metaclass=ABCMeta):
    def __init__(self, timeout: int = 10):
        super().__init__(timeout)
        self.headers = None
        self.read_count = 0
        self._writer = None

    def _finish_thread(self) -> None:
        self._writer.flush()
        self._writer.close()

    def _write(self, data: bytes) -> None:
        self._writer.write(data)

    def download_headers(
        self, url: Union[str, ParseResult], *ticket_args, **ticket_kwargs
    ) -> dict:
        """
        Downloads and returns only the headers.

        The only reliable way to get the headers without downloading the entire
        BAM file is to use the new v1.2 'class' keyword, which annotates a URL as
        containing 'header'.

        Args:
            url:
            *ticket_args:
            **ticket_kwargs:

        Returns:
            The header dict. Has the side-effect of initializing self.headers.
        """
        if "class" not in ticket_kwargs:
            ticket_kwargs["class"] = HEADER_CLASS

        ticket = get_ticket(url, *ticket_args, **ticket_kwargs)

        # If 'header' URLs are provided, they must come before 'body' URLs.
        if "class" in ticket["urls"][0]:
            header_urls = [
                url_object
                for url_object in ticket["urls"]
                if url_object["class"] == HEADER_CLASS
            ]
            return self.download_header_urls(header_urls)
        else:
            raise DownloaderStateError(
                "Headers have not yet been initialized, and the server does not "
                "provide header-only URLs."
            )

    def download_header_urls(self, header_url_objects) -> dict:
        """
        Downloads only the headers from the given url objects.

        Threading is not necessary for this - we simply download (and concatenate if
        necessary) all contents of the URLs into a byte string, create a
        pysam.AlignmentFile, and extract the headers.

        Args:
            header_url_objects:

        Returns:
            The header dict. Has the side-effect of initializing self.headers.
        """
        bam_bytes = io.BytesIO()

        for url_object in header_url_objects:
            handle_htsget_url(
                parse_url(url_object[URL_KEY]),
                url_object,
                self.timeout,
                bam_bytes.write,
            )

        pipe_reader, pipe_writer = os.pipe()

        with open(pipe_writer, "wb") as out:
            out.write(bam_bytes.getvalue())
            out.write(BGZF_EOF_MARKER)
            out.flush()

        with open(pipe_reader, "rb") as inp:
            bam = pysam.AlignmentFile(inp, "r", check_sq=False)
            self.headers = bam.header
            return self.headers

    def __iter__(self) -> Iterator[SamRecord]:
        pass


class BamHtsgetDownloader(HtsgetDownloader):
    """
    Downloader that only saves BAM data to a file.
    """
    def __init__(self, output_file: Path, timeout: int = 10):
        super().__init__(timeout)
        self.byte_count = 0
        self.output_file = output_file
        self._output_fh = None

    def _init_thread(self):
        self._output_fh = open(self.output_file, "wb")

    def _finish_thread(self) -> None:
        self._output_fh.flush()
        self._output_fh.close()

    def _write(self, data: bytes) -> None:
        self._output_fh.write(data)
        self.byte_count += len(data)


def get_ticket(
    url: Union[str, ParseResult], *ticket_args, timeout: int = 10, **ticket_kwargs
) -> dict:
    if isinstance(url, str):
        parsed_url = parse_url(url)
    else:
        parsed_url = cast(ParseResult, url)

    ticket_request_url = get_ticket_request_url(
        parsed_url, *ticket_args, **ticket_kwargs
    )
    logging.debug(f"handle_ticket_request(url={ticket_request_url})")
    response = httpget(ticket_request_url, timeout=timeout)
    response_json = response.json()

    if HTSGET_KEY not in response_json:
        raise ValueError(f"Invalid response: {response_json}")

    return response_json[HTSGET_KEY]


def get_ticket_request_url(
    parsed_url: ParseResult,
    data_format: Optional[str] = None,
    class_: Optional[str] = None,
    reference_name: Optional[str] = None,
    reference_md5: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    # It's not clear to me what happens if I request BAM data with only a subset
    # of fields. Excluding this for now.
    # fields: Sequence[str] = ('QNAME', 'FLAGS', 'SEQ', 'QUAL'),
    tags: Optional[Sequence[str]] = None,
    notags: Optional[Sequence[str]] = None,
    **kwargs
) -> str:
    """
    Generates an htsget request URL.

    Args:
        parsed_url: The URL of the data to retrieve, parsed by urlparse.
        data_format: The requested format of the returned data.
        class_: Class of data to request; by default, complete data is queried.
            Currently, the only valid value to this parameter is 'header'.
        reference_name: The reference sequence name, for example "chr1",
            "1", or "chrX". If unspecified, all data is returned.
        reference_md5: The MD5 checksum uniquely representing the reference
            sequence as a lower-case hexadecimal string, calculated as the MD5
            of the upper-case sequence excluding all whitespace characters
            (this is equivalent to SQ:M5 in SAM).
        start: The start position of the range on the reference, 0-based,
            inclusive. If specified, ``reference_name`` or ``reference_md5``
            must also be specified.
        end: The end position of the range on the reference, 0-based
            exclusive. If specified, ``reference_name`` or ``reference_md5``
            must also be specified.
        tags: Sequence of tags to include.
        notags: Seqeunce of tags to exclude

    Todo:
        error checking

    Returns:
        A URL.
    """
    get_vars: Dict[Any, Union[Any, List[Any]]] = parse_qs(parsed_url.query)

    if data_format is not None:
        get_vars["format"] = data_format.upper()

    if not class_ and "class" in kwargs:
        class_ = kwargs["class"]

    if class_ is not None:
        if class_ != HEADER_CLASS:
            raise ValueError("'class' must be 'header'")
        get_vars["class"] = class_

    if reference_name is not None:
        get_vars["referenceName"] = reference_name

    if reference_md5 is not None:
        get_vars["referenceMD5"] = reference_md5

    if start is not None:
        get_vars["start"] = int(start)

    if end is not None:
        get_vars["end"] = int(end)

    # if fields is not None:
    #     get_vars["fields"] = ",".join(fields)

    if tags is not None:
        get_vars["tags"] = ",".join(tags)

    if notags is not None:
        get_vars["notags"] = ",".join(notags)

    new_url = list(parsed_url)
    new_url[4] = urlencode(get_vars, doseq=True)

    return urlunparse(new_url)


def handle_htsget_url(
    parsed_url: ParseResult,
    url_object: dict,
    timeout: int,
    writer: Callable[[bytes], None],
):
    if parsed_url.scheme.startswith("http"):
        _handle_http_url(
            url_object[URL_KEY], url_object.get("headers", ""), timeout, writer
        )
    elif parsed_url.scheme == "data":
        _handle_data_uri(parsed_url, writer)
    else:
        raise ValueError("Unsupported URL scheme: {}".format(parsed_url.scheme))


def _handle_http_url(
    url: str, headers: str, timeout: int, writer: Callable[[bytes], None]
) -> str:
    logging.debug(f"handle_http_url(url={url}, headers={headers})")
    response = httpget(url, headers=headers, stream=True, timeout=timeout)
    length = 0
    piece_size = DEFAULT_CHUNK_SIZE

    for piece in response.iter_content(piece_size):
        length += len(piece)
        writer(piece)

    if CONTENT_LENGTH in response.headers:
        content_length = int(response.headers[CONTENT_LENGTH])

        if content_length != length:
            raise ContentLengthMismatchError(
                "Length mismatch {content_length} != {length}"
            )

    # Assume lowest-common-denominator
    htsget_version = "1.0"

    if CONTENT_TYPE in response.headers:
        content_type = response.headers[CONTENT_TYPE]
        content_type_match = CONTENT_TYPE_RE.match(content_type)
        if content_type_match:
            htsget_version = content_type_match.group(1)

    return htsget_version


def _handle_data_uri(parsed_url: ParseResult, writer: Callable[[bytes], None]):
    split = parsed_url.path.split(",", 1)
    # TODO parse out the encoding properly.
    description = split[0]
    data = base64.b64decode(split[1])
    logging.debug(f"handle_data_uri({description}, length={len(data)})")
    writer(data)


def httpget(*args, **kwargs):
    """
    Send a GET request and return the response.
    """
    response = requests.get(*args, **kwargs)
    response.raise_for_status()
    return response
