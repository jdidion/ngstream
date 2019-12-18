import pkg_resources

from ngstream import api
from ngstream.api import Protocol, Record, Fragment, dump_fastq
from ngstream.utils import GenomeReference


try:
    __version__ = pkg_resources.get_distribution(__name__).version
except pkg_resources.DistributionNotFound:
    __version__ = "Unknown"


def open(accession: str, protocol: str, start: bool = True, **kwargs) -> Protocol:
    """
    Opens a streaming protocol.

    Args:
        accession: The accession of the dataset; this is protocol-specific - for 'sra'
            it is a string identifier, and for htsget it is a URL string.
        protocol: Name of the data transfer protocol. Currently 'sra' and 'htsget'
            are supported.
        start: Whether to call `protocol.start()`.
        kwargs: Additional arguments to pass to the :ngstream.api.Protocol: class;
            see the protocol documentation for details.

    Returns:
        A :ngstream.api.Protocol: object, which is iterable and a context manager.
    """
    protocol_cls = api.get_protocol(protocol)
    protocol = protocol_cls(accession, **kwargs)

    if start:
        protocol.start()

    return protocol
