# -*- coding: utf-8 -*-
"""
"""
from types import FileArg, BinMode
from xphyle.utils import fileinput

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
        self.file_format = file_format
    
    @property
    def delivers_qualities(self):
        return self.file_format.delivers_qualities

class FileSeqIO(FormatSeqIO):
    """Base class for SeqIO classes that read from a file.

    Args:
        files: Path or file-like object. The file may be compressed
            (.gz, .bz2, .xz)
        mode: The file open mode. Must be binary.
        file_format: A file format name, or an instance of SeqFileFormat
        kwargs: Additional arguments to pass to open_
    """
    def __init__(self, *files: FileArg, mode: str = 'b',
                 file_format: SequenceFormat, **kwargs):
        if 'b' not in mode:
            raise ValueError("'mode' must be binary")
        super(FileSeqIO, self).__init__(file_format)
        self.reader = self._open_reader(*files, mode=mode, **kwargs)
    
    def _open_reader(self, *files: FileArg, mode: str = 'b', **kwargs):
        return fileinput(files, BinMode)
    
    @property
    def name(self):
        return self.fileobj.name
    
    def close(self):
        self.fileobj.close()

class SingleReader(object):
    paired = False

class SingleFileReader(FileSeqIO, SingleReader):
    def __init__(self, *files, file_format, **kwargs):
        super(SingleFileReader, self).__init__(
            *files, 'rb', file_format, **kwargs)
    
    def __next__(self):
        return self.file_format.read_record(self.reader)

class PairedReader(object):
    paired = True
    
    def create_record(self, reads):
        if reads[0].name != reads[1].name:
            raise FormatError(
                "Reads in pair do not have same name: ({0!r} != {1!r})".format(
                reads[0].name, reads[1].name))
        return reads

class PairedFileReader(FormatSeqIO, PairedReader):
    """Read from a pair of (possibly compressed) files containing sequences.

    Args:
        name: A name for this sequence reader
        read1, read2: SeqIO instances
        file_format: A file format name, or an instance of SeqFileFormat
    """
    def __init__(self, name, read1, read2, file_format):
        super(PairedFileReader, self).__init__(file_format)
        self.name = name
        self.read1 = read1
        self.read2 = read2
    
    def __next__(self):
        reads = tuple(next(self.read1), next(self.read2))
        if not all(reads):
            raise FormatError("At least one read is 'None'")
        return self.create_record(reads)
    
    def close(self):
        self.read1.close()
        self.read2.close()

class InterleavedFileReader(FileSeqIO, PairedReader):
    def __init__(self, files: FileArg, file_format, **kwargs):
        super(InterleavedFileReader, self).__init__(
            path, 'rb', file_format, **kwargs)
    
    def __next__(self):
        return self.create_record(self.file_format.read_pair(self.fileobj))
    
    def iter_single_end(self, end):
        end -= 1
        for reads in itr(self):
            yield reads[end]

class SingleWriter(object):
    paired = False

class SequenceWriter(FileSeqIO):
    def __init__(self, path, file_format, **kwargs):
        super(SequenceWriter, self).__init__(
            path, 'wb', file_format, **kwargs)
    
    def write(self, record):
        self.fileobj.write(self.file_format.format_record(record))

class PairedFileWriter(FormatSeqIO):
    """Write sequences to a pair of (possibly compressed) files.

    Args:
        name: A name for this sequence reader
        read1, read2: SeqIO instances
        file_format: A file format name, or an instance of SeqFileFormat
    """
    paired = True
    
    def __init__(self, name, read1, read2, file_format):
        super(PairedFileWriter, self).__init__(file_format)
        self.name = name
        self.read1 = read1
        self.read2 = read2
    
    def write(self, read1, read2):
        self.read1.write(read1)
        self.read2.write(read2)
    
    def close(self):
        self.read1.close()
        self.read2.close()

class InterleavedFileWriter(FileSeqIO):
    paired = True
    
    def __init__(self, path, file_format, **kwargs):
        super(InterleavedFileWriter, self).__init__(
            path, 'wb', file_format, **kwargs)
    
    def write(self, read1, read2):
        self.fileobj.write(self.file_format.format_record(read1))
        self.fileobj.write(self.file_format.format_record(read2))
