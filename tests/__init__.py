import os
import urllib

def no_internet():
    """Test whether there's no internet connection available.
    """
    try:
        urllib.request.urlopen("https://github.com").info()
        return False
    except:
        return True

def load_expected(fname):
    path = os.path.join(os.path.dirname(__file__), 'expected', fname)
    with open(path, 'rt') as inp:
        return inp.read()