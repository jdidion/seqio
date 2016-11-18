# -*- coding: utf-8 -*-
"""
"""
from seqio.formats import get_file_format

# Exceptions

class FormatError(Exception):
    """Raised when an input file (FASTA or FASTQ) is malformatted.
    """
    pass

# Base classes

class SeqIO(object):
    """Base class for classes that read/write sequences.
    
    Subclasses must provide:
    * A member 'name'
    * A function 'close(self)'
    """
    def __enter__(self):
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        self.close()
    
    def __repr__(self):
        return "<{0!r}(name={1!r})>".format(self.class, self.name)

class FormatSeqIO(SeqIO):
    """Base class for SeqIO classes with a specific file format.
    
    Args:
        file_format: The file format
    """
    def __init__(self, file_format):
        if isinstance(file_format, str):
            file_format = get_file_format(file_format)
        self.file_format = file_format
    
    @property
    def delivers_qualities(self):
        return self.file_format.delivers_qualities

class FileSeqIO(FormatSeqIO):
    """Base class for SeqIO classes that read from a file.

    Args:
        path: Path or file-like object. The file may be compressed
            (.gz, .bz2, .xz)
        mode: The file open mode. Must be binary.
        file_format: A file format name, or an instance of SeqFileFormat
        kwargs: Additional arguments to pass to open_
    """
    def __init__(self, path, mode, file_format, **kwargs):
        if 'b' not in mode:
            raise ValueError("'mode' must be binary")
        super(FileSeqIO, self).__init__(file_format)
        self.fileobj = self.file_format.open(path, mode=mode, **kwargs)
    
    @property
    def name(self):
        return self.fileobj.name
    
    def close(self):
        self.fileobj.close()
