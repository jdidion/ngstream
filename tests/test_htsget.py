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

from ngstream.protocols.htsget import HtsgetProtocol
from ngstream.writers import BufferWriter
from ngstream.formats import FastqFileFormat
from ngstream.utils import CoordinateBatcher, GenomeReference

PORT = 6160
SERVER_URL = "http://localhost:{}".format(PORT)
REFERENCE = dict(chr1=100)


class MockURLInstance:
    root = os.path.dirname(__file__)

    def __init__(
        self,
        input_file,
        expected_file1,
        expected_file2=None,
        headers={},
        error_code=None,
        truncate=False,
    ):
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
            with open(os.path.join(self.root, "input", self.input_file), "rb") as inp:
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
        self.reference = GenomeReference("mock", REFERENCE)


class TestDataTransfers(ServerTest):
    """
    Test cases for various data transfers.
    """

    def assert_data_transfer_ok(self, test_instances, max_retries=0):
        self.httpd.test_instances = test_instances
        reader = HtsgetReader(MockRequestHandler.ticket_url, self.reference)
        buf = BufferWriter(True)
        with reader, FastqFileFormat(buf) as writer:
            for read in reader:
                print(read)
                writer(*read)
        all_data = [
            "".join(t.expected[i] for t in test_instances if t.expected[i] is not None)
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
        instances = [MockURLInstance("paired.bam", "paired.1.fastq", "paired.2.fastq")]
        self.assert_data_transfer_ok(instances)


import os.path
import argparse
import random
import sys
import logging
import multiprocessing
import traceback

import tqdm
import hjson
import humanize
import pandas as pd

import tester


class TestCase(object):
    """
    A single client/server combination for a data file.
    """
    def __init__(self, input_file, server_name, url, client):
        self.input_file = input_file
        self.server_name = server_name
        self.url = url
        self.client = client
        self.random_seed = 1
        self.max_references = 100
        self.num_random_queries = 5
        self.tmpdir = None
        self.filter_unmapped = False
        self.log_file_prefix = None

    @property
    def log_file_name(self):
        return os.path.join(self.log_file_prefix, "{}.log".format(self.test_name))

    @property
    def test_name(self):
        data_name = os.path.split(self.input_file)[-1]
        return "{}_{}_{}".format(data_name, self.server_name, self.client)


def parse_input(json_cases, data_file_prefix, clients):
    """
    Parses the input JSON into test cases.
    """
    cases = []
    for json_case in json_cases["cases"]:
        input_file = os.path.join(data_file_prefix, json_case["file"])
        for json_server in json_case["servers"]:
            for client in clients:
                cases.append(TestCase(
                    input_file, json_server["name"], json_server["url"], client))
                if "filter_unmapped" in json_server:
                    cases[-1].filter_unmapped = json_server["filter_unmapped"]
    return cases


def run_test(case):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    log_handler = logging.FileHandler(case.log_file_name, mode="w")
    formatter = logging.Formatter('%(asctime)s %(message)s')
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    random.seed(case.random_seed)
    server_tester = tester.ServerTester(
        case.input_file, case.url, filter_unmapped=case.filter_unmapped,
        tmpdir=case.tmpdir, max_references=case.max_references,
        client=tester.client_map[case.client])
    completed = False
    try:
        server_tester.initialise()
        server_tester.run_full_contig_fetch()
        server_tester.run_start_reads()
        server_tester.run_end_reads()
        server_tester.run_random_reads(case.num_random_queries)
        logger.info("Tests completed")
        completed = True
    except tester.DownloadFailedException as dfe:
        print("Download failed: ", dfe, "see ", case.log_file_name, file=sys.stderr)
        logger.info("Download failed: {}".format(dfe))
    except tester.TestFailedException as tfe:
        print(
            "Test failed:", tfe, ". Is the input file correct?", case.log_file_name,
            file=sys.stderr)
        logger.info("Test failed: {}".format(tfe))
    except Exception as e:
        logger.info("Unexpected error on {}: {}".format(case.log_file_name, e))
        print("Unexpected error! See", case.log_file_name)
        traceback.print_exc()
    finally:
        server_tester.cleanup()
        logger.removeHandler(log_handler)
        log_handler.close()
    ret = None
    if completed:
        ret = server_tester.report()
    return case, ret
