import pytest
import urllib


try:
    urllib.request.urlopen("https://github.com").info()
    INTERNET_CONNECTION = True
except:
    INTERNET_CONNECTION = False


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "skipif_no_internet: mark test to run only if there is an internet connection")


def pytest_runtest_setup(_):
    if not INTERNET_CONNECTION:
        pytest.skip("test requires an internet connection")
