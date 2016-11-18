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





class FastaQualReader(object):
    """
    Reader for reads that are stored in .(CS)FASTA and .QUAL files.
    """
    delivers_qualities = True

    def __init__(self, fastafile, qualfile, sequence_class=Sequence):
        """
        fastafile and qualfile are filenames or file-like objects.
        If a filename is used, then .gz files are recognized.

        The objects returned when iteritng over this file are instances of the
        given sequence_class.
        """
        self.fastareader = FastaReader(fastafile)
        self.qualreader = FastaReader(qualfile, keep_linebreaks=True)
        self.sequence_class = sequence_class

    def __iter__(self):
        """
        Yield Sequence objects.
        """
        # conversion dictionary: maps strings to the appropriate ASCII-encoded character
        conv = dict()
        for i in range(-5, 256 - 33):
            conv[str(i)] = chr(i + 33)
        for fastaread, qualread in zip(self.fastareader, self.qualreader):
            if fastaread.name != qualread.name:
                raise FormatError("The read names in the FASTA and QUAL file "
                    "do not match ({0!r} != {1!r})".format(fastaread.name, qualread.name))
            try:
                qualities = ''.join([conv[value] for value in qualread.sequence.split()])
            except KeyError as e:
                raise FormatError("Within read named {0!r}: Found invalid quality "
                    "value {1}".format(fastaread.name, e))
            assert fastaread.name == qualread.name
            yield self.sequence_class(fastaread.name, fastaread.sequence, qualities)

    def close(self):
        self.fastareader.close()
        self.qualreader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
