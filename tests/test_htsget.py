"""
Test cases using a simple HTTP server running locally.
Lots of code borrowed from Jerome Kelleher's htsget python library, which is
under the Apache 2.0 license.
"""
from unittest import TestCase
import json
import os
import sys
import tempfile
import threading
from unittest import mock
from http.server import BaseHTTPRequestHandler
import socketserver
from urllib.parse import urljoin, urlparse

from ngstream.htsget import HtsgetReader
from ngstream.writers import FastqWriter, BufferWriter
from ngstream.utils import CoordinateBatcher, GenomeReference

PORT = 6160
SERVER_URL = "http://localhost:{}".format(PORT)
REFERENCE = dict(chr1=100)

class MockURLInstance():
    root = os.path.dirname(__file__)
    
    def __init__(
            self, input_file, expected_file1, expected_file2=None, headers={}, 
            error_code=None, truncate=False):
        self.url = "/input/{}".format(input_file)
        self.input_file = input_file
        self.expected_files = [
            os.path.join(self.root, "expected", fname)
            for fname in (expected_file1, expected_file2)
            if fname is not None
        ]
        self._data = None
        self._expected = [None, None]
        self.headers = headers
        self.error_code = error_code
        self.truncate = truncate
    
    @property
    def data(self):
        if self._data is None:
            with open(os.path.join(self.root, 'input', self.input_file), 'rb') as inp:
                self._data = inp.read()
        return self._data
    
    @property
    def expected(self):
        if all(d is None for d in self._data):
            for i, fname in enumerate(self.expected):
                if fname is not None:
                    with open(fname, 'rt') as inp:
                        self._data[i] = inp.read()
        return self._data


class MockServer(socketserver.TCPServer):
    """
    A local test server designed to be run in a thread and shutdown smoothly
    for test cases.
    """
    allow_reuse_address = True
    # This is set by clients to contain the list of MockURLInstance objects.
    test_instances = []

    def shutdown(self):
        self.socket.close()
        socketserver.TCPServer.shutdown(self)


class MockRequestHandler(BaseHTTPRequestHandler):
    ticket_path = "/ticket"
    ticket_url = urljoin(SERVER_URL, ticket_path)

    def log_message(self, format, *args):
        # Silence the logger.
        pass

    def do_GET(self):
        url_map = {instance.url: instance for instance in self.server.test_instances}
        parsed = urlparse(self.path)
        if parsed.path == self.ticket_path:
            self.send_response(200)
            self.end_headers()
            urls = [
                {
                    "url": urljoin(SERVER_URL, test_instance.url),
                    "headers": test_instance.headers
                } for test_instance in self.server.test_instances]
            ticket = {
                "urls": urls
            }
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


class ServerTest(TestCase):
    """
    Superclass of tests needing a server running in a thread.
    """
    @classmethod
    def setup_class(cls):
        cls.httpd = MockServer(("", PORT), MockRequestHandler)

        def target():
            # Sometimes the server doesn't shutdown cleanly, but we don't
            # care here.
            try:
                cls.httpd.serve_forever()
            except ValueError:
                pass

        cls.httpd_thread = threading.Thread(target=target)
        cls.httpd_thread.start()

    @classmethod
    def teardown_class(cls):
        cls.httpd.shutdown()
        cls.httpd_thread.join()

    def setUp(self):
        self.reference = GenomeReference('mock', REFERENCE)


class TestDataTransfers(ServerTest):
    """
    Test cases for various data transfers.
    """
    def assert_data_transfer_ok(self, test_instances, max_retries=0):
        self.httpd.test_instances = test_instances
        reader = HtsgetReader(
            MockRequestHandler.ticket_url, self.reference)
        buf = BufferWriter(True)
        with reader, FastqWriter(buf) as writer:
            for read in reader:
                print(read)
                writer(*read)
        all_data = [
            "".join(
                t.expected[i] 
                for t in test_instances
                if t.expected[i] is not None)
            for i in (0, 1)
        ]
        self.assertEqual(buf.value, all_data)

    # def test_simple_data(self):
    #     instances = [
    #         MockURLInstance(url="/data1", data=b"data1"),
    #         MockURLInstance(url="/data2", data=b"data2")
    #     ]
    #     self.assert_data_transfer_ok(instances)

    def test_binary_data(self):
        instances = [
            MockURLInstance('paired.bam', 'paired.1.fastq', 'paired.2.fastq')
        ]
        self.assert_data_transfer_ok(instances)
