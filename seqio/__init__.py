# -*- coding: utf-8 -*-
"""The main seqio methods.
"""
from seqio.sequence import *

class FormatError(Exception):
    """Raised when an input file (FASTA or FASTQ) is malformatted.
    """
    pass

class SeqIO(object):
    """Base class for classes that read/write sequences.
    
    Subclasses must provide:
    * A member 'name'
    * A function 'close(self)'
    """
    def __init__(self, file_format):
        if isinstance(file_format, str):
            file_format = get_file_format(file_format)
        self.file_format = file_format
    
    @property
    def delivers_qualities(self):
        return self.file_format.delivers_qualities
        
    def __enter__(self):
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        self.close()
    
    def __repr__(self):
        return "<{0!r}(name={1!r})>".format(self.class, self.name)

class FileSeqIO(SeqIO):
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

class SequenceReader(FileSeqIO):
    paired = False
    
    def __init__(self, path, file_format, **kwargs):
        super(FileSequenceReader, self).__init__(
            path, 'rb', file_format, **kwargs)
    
    def __next__(self):
        return self.file_format.read_record(self.fileobj)

class PairedSeqIO(SeqIO):
    """Open a pair of possibly compressed file containing sequences.

    Args:
        name: A name for this sequence reader
        read1, read2: SeqIO instances
        file_format: A file format name, or an instance of SeqFileFormat
    """
    def __init__(self, name, read1, read2, file_format):
        super(PairedSeqIO, self).__init__(file_format)
        self.name = name
        self.read1 = read1
        self.read2 = read2
    
    def close(self):
        self.read1.close()
        self.read2.close()

class PairedReader(PairedSeqIO):
    paired = True
    
    def __next__(self):
        read1 = next(self.read1)
        read2 = next(self.read2)
        if read1 and read2:
            return self.make_record(read1, read2)
        else:
            raise FormatError("At least one read is 'None'")
    
    def make_record(self, read1, read2):
        return (read1, read2)

class InterleavedReader(SequenceReader):
    paired = True
    
    def __next__(self):
        return self.file_format.read_pair(self.fileobj)
    
    def read_single_end(self, end):
        self.paired = False
        if end == 1:
            self.make_record = lambda self: \
                self.file_format.read_pair(self.fileobj)[0]
        else:
            self.make_record = lambda self: \
                self.file_format.read_pair(self.fileobj)[1]




class SequenceWriter(object):
    @property
    def paired(self):
        raise NotImplemented()
    
    def write(self, record):
        raise NotImplemented()
