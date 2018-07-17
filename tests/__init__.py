import tempfile


def mktmpfile(tmpdir, **kwargs):
    return tempfile.mkstemp(dir=tmpdir, **kwargs)[1]
