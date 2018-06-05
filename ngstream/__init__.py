"""ngstream user interface.
"""
from ngstream import api
from ngstream.api import Protocol, Record, Fragment
from ngstream.utils import GenomeReference
from ngstream._version import get_versions
__version__ = get_versions()['version']
del get_versions


def protocol(name: str, accn: str, **kwargs) -> Protocol:
    """Open a streaming protocol.

    Args:
        name: Name of the data transfer protocol. Currently 'sra' and 'htsget'
            are supported.
        accn: The accession of the dataset; this is protocol-specific - for 'sra' it is
            a string identifier, and for htsget it is a URL string.
        **kwargs: Additional arguments to pass to the :ngstream.api.Protocol: class;
            see the protocol documentation for details.

    Returns:
        A :ngstream.api.Protocol: object, which is iterable and a context manager.
    """
    protocol_cls = api.get_protocol(name)
    return protocol_cls(accn, **kwargs)
