# -*- coding: utf-8 -*-
"""
"""
from seqio.base import *
from seqio.sequence import *

class SingleFileReader(FileSeqIO):
    paired = False
    
    def __init__(self, path, file_format, **kwargs):
        super(SequenceReader, self).__init__(
            path, 'rb', file_format, **kwargs)
    
    def __next__(self):
        return self.file_format.read_record(self.fileobj)

class PairedReader(object):
    paired = True
    
    def create_record(self, reads):
        if reads[0].name != reads[1].name:
            raise FormatError(
                "Reads in pair do not have same name: ({0!r} != {1!r})".format(
                reads[0].name, reads[1].name))
        return reads

class PairedFileReader(SeqIO, PairedReader):
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
    def __init__(self, path, file_format, **kwargs):
        super(InterleavedFileReader, self).__init__(
            path, 'rb', file_format, **kwargs)
    
    def __next__(self):
        return self.create_record(self.file_format.read_pair(self.fileobj))
    
    def iter_single_end(self, end):
        end -= 1
        for reads in itr(self):
            yield reads[end]
