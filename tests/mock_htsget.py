from http.server import BaseHTTPRequestHandler
import json
import os
from pathlib import Path
import socketserver
from typing import Sequence
from urllib.parse import urljoin, urlparse


PORT = 6160
SERVER_URL = "http://localhost:{}".format(PORT)
TICKET_PATH = "/ticket"
TICKET_URL = urljoin(SERVER_URL, TICKET_PATH)
REFERENCE = dict(chr1=100)


class MockURLInstance:
    def __init__(
        self,
        root: Path,
        input_file,
        expected_file1,
        expected_file2=None,
        headers=None,
        error_code=None,
        truncate=False,
    ):
        self.root = root
        self.url = f"/input/{input_file}"
        self.input_file = self.root / input_file
        self.expected_files = [
            os.path.join(self.root, "expected", fname)
            for fname in (expected_file1, expected_file2)
            if fname is not None
        ]
        self._data = None
        self._expected = [None, None]
        self.headers = headers or {}
        self.error_code = error_code
        self.truncate = truncate

    @property
    def data(self):
        if self._data is None:
            with open(self.root / "input" / self.input_file, "rb") as inp:
                self._data = inp.read()
        return self._data

    @property
    def expected(self):
        if all(d is None for d in self._data):
            for i, fname in enumerate(self.expected):
                if fname is not None:
                    with open(fname, "rt") as inp:
                        self._data[i] = inp.read()
        return self._data


class MockServer(socketserver.TCPServer):
    """
    A local test server designed to be run in a thread and shutdown smoothly
    for test cases.
    """
    allow_reuse_address = True
    # This is set by clients to contain the list of MockURLInstance objects.
    test_instances: Sequence[MockURLInstance] = []

    def shutdown(self):
        self.socket.close()
        socketserver.TCPServer.shutdown(self)


class MockRequestHandler(BaseHTTPRequestHandler):
    def log_message(self, _format, *args):
        # Silence the logger.
        pass

    def do_GET(self):
        url_map = {instance.url: instance for instance in self.server.test_instances}
        parsed = urlparse(self.path)
        if parsed.path == TICKET_PATH:
            self.send_response(200)
            self.end_headers()
            urls = [
                {
                    "url": urljoin(SERVER_URL, test_instance.url),
                    "headers": test_instance.headers,
                }
                for test_instance in self.server.test_instances
            ]
            ticket = {"urls": urls}
            self.wfile.write(json.dumps(ticket).encode())
        elif self.path in url_map:
            instance = url_map[self.path]
            if instance.error_code is not None:
                self.send_error(instance.error_code)
            else:
                self.send_response(200)
                self.send_header("Content-Length", len(instance.data))
                if instance.truncate:
                    self.end_headers()
                    self.wfile.write(instance.data[:-1])
                else:
                    self.end_headers()
                    self.wfile.write(instance.data)
        else:
            self.send_error(404)
