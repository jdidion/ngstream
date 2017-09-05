"""
Test cases using a simple HTTP server running locally.
Lots of code borrowed from Jerome Kelleher's htsget python library, which is
under the Apache 2.0 license.
"""
import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest import mock
from http.server import BaseHTTPRequestHandler
import time
import socketserver
from urllib.parse import urljoin

from ngstream.htsget import HtsgetReader
from ngstream.writers import FastqWriter, BufferWriter
from ngstream.utils import CoordinateBatcher, GenomeReference

PORT = 6160
SERVER_URL = "http://localhost:{}".format(PORT)
REFERENCE = dict(chr1=100)

class TestUrlInstance(object):
    def __init__(self, url, data, headers={}, error_code=None, truncate=False):
        self.url = url
        self.data = data
        self.headers = headers
        self.error_code = error_code
        self.truncate = truncate


class TestServer(socketserver.TCPServer):
    """
    A local test server designed to be run in a thread and shutdown smoothly
    for test cases.
    """
    allow_reuse_address = True
    # This is set by clients to contain the list of TestUrlInstance objects.
    test_instances = []

    def shutdown(self):
        self.socket.close()
        socketserver.TCPServer.shutdown(self)


class TestRequestHandler(BaseHTTPRequestHandler):

    ticket_path = "/ticket"
    ticket_url = urljoin(SERVER_URL, ticket_path)

    def log_message(self, format, *args):
        # Silence the logger.
        pass

    def do_GET(self):
        url_map = {instance.url: instance for instance in self.server.test_instances}
        if self.path == self.ticket_path:
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


class ServerTest(unittest.TestCase):
    """
    Superclass of tests needing a server running in a thread.
    """
    @classmethod
    def setup_class(cls):
        cls.httpd = TestServer(("", PORT), TestRequestHandler)

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
            TestRequestHandler.ticket_url, self.reference)
        buf = BufferWriter(True)
        with reader, FastqWriter(buf) as writer:
            for read in reader:
                writer(*read)
        all_data = b"".join(test_instance.data for test_instance in test_instances)
        self.assertEqual(buf.value, all_data)

    def test_simple_data(self):
        instances = [
            TestUrlInstance(url="/data1", data=b"data1"),
            TestUrlInstance(url="/data2", data=b"data2")
        ]
        self.assert_data_transfer_ok(instances)

    def test_binary_data(self):
        instances = []
        for j in range(10):
            instances.append(TestUrlInstance(
                url="/path/to/data/{}".format(j),
                data=bytes(j) * 1024))
        self.assert_data_transfer_ok(instances)

    def test_transfer_with_cli(self):
        test_instances = [
            TestUrlInstance(url="/data1", data=b"data1"),
            TestUrlInstance(url="/data2", data=b"data2")
        ]
        self.httpd.test_instances = test_instances
        try:
            fd, filename = tempfile.mkstemp()
            os.close(fd)
            cmd = [TestRequestHandler.ticket_url, "-O", filename]
            parser = cli.get_htsget_parser()
            args = parser.parse_args(cmd)
            with mock.patch("sys.exit") as mocked_exit:
                cli.run(args)
                mocked_exit.assert_called_once_with(0)
            all_data = b"".join(test_instance.data for test_instance in test_instances)
            with open(filename, "rb") as f:
                self.assertEqual(f.read(), all_data)
        finally:
            os.unlink(filename)

    def test_transfer_with_cli_stdout(self):
        test_instances = [
            TestUrlInstance(url="/data1", data=b"data1"),
            TestUrlInstance(url="/data2", data=b"data2")
        ]
        self.httpd.test_instances = test_instances
        saved = sys.stdout
        sys.stdout = self.output_file
        try:
            cmd = [TestRequestHandler.ticket_url]
            parser = cli.get_htsget_parser()
            args = parser.parse_args(cmd)
            with mock.patch("sys.exit") as mocked_exit:
                cli.run(args)
                mocked_exit.assert_called_once_with(0)
            all_data = b"".join(test_instance.data for test_instance in test_instances)
            self.output_file.seek(0)
            self.assertEqual(self.output_file.read(), all_data)
        finally:
            sys.stdout = saved
