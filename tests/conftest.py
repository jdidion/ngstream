import json
from pathlib import Path
import pytest
import threading
from typing import Sequence, cast
from urllib.request import urlopen
from ngstream.utils import GenomeReference
from .mock_htsget import MockServer, MockRequestHandler, PORT, REFERENCE
from .real_htsget import ServerTester

try:
    urlopen("https://github.com").info()
    INTERNET_CONNECTION = True
except:
    INTERNET_CONNECTION = False


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "skipif_no_internet: mark test to run only if there is an internet connection")


def pytest_runtest_setup(item):
    if list(item.iter_markers(name='skipif_no_internet')) and not INTERNET_CONNECTION:
        pytest.skip("test requires an internet connection")


@pytest.fixture(scope='session')
def server():
    httpd = MockServer(("", PORT), MockRequestHandler)

    def target():
        # Sometimes the server doesn't shutdown cleanly, but we don't
        # care here.
        try:
            httpd.serve_forever()
        except ValueError:
            pass

    httpd_thread = threading.Thread(target=target)
    httpd_thread.start()

    try:
        yield httpd
    finally:
        httpd.shutdown()
        httpd_thread.join()


@pytest.fixture
def reference():
    yield GenomeReference("mock", REFERENCE)


@pytest.fixture
def htsget_json(datadir: Path) -> Sequence[dict]:
    with open(datadir / 'servers.json', 'rt') as inp:
        return json.load(inp)['cases']


@pytest.fixture(params=htsget_json)
def htsget_case(request) -> ServerTester:
    for case in cast(Sequence['dict'], request.param):
        filename = case['file']
        for _server in case['servers']:
            yield ServerTester(
                filename, _server['name'], _server['url'],
                _server.get('filter_unmapped', False),
                int(_server.get('random_seed', 1)),
                _server.get('bearer_token', None),
                _server.get('ca_bundle', None))
