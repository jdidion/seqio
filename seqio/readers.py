# -*- coding: utf-8 -*-
"""
"""
from seqio.base import *
from seqio.sequence import *

class SingleReader(object):
    paired = False

class SingleFileReader(FileSeqIO, SingleReader):
    def __init__(self, path, file_format, **kwargs):
        super(SequenceReader, self).__init__(
            path, 'rb', file_format, **kwargs)
    
    def __next__(self):
        return self.file_format.read_record(self.fileobj)

class FastaQualReader(SeqIO, SingleReader):
    """Reader for reads that are stored in .(CS)FASTA and .QUAL files.
    
    Args:
        fastafile: path or file-like object of sequences in FASTA format
        qualfile: path or file-like object of qualities in FASTA format
        fasta_format: An instance of ``seqio.formats.Fasta``
        kwargs: Additional arguments passed to the file open method
    """
    delivers_qualities = True
    
    def __init__(self, fastafile, qualfile, sequence_class=Sequence, **kwargs):
        self.fasta_reader = SingleFileReader(fastafile, Fasta(), **kwargs)
        self.qual_reader = FastaReader(qualfile, Fasta(linesep=b' '), **kwargs)
        self.sequence_class = sequence_class
        self.qual_conversion = {}
        for i in range(-5, 256 - 33):
            self.qual_conversion[str(i).encode()] = chr(i + 33).encode()
    
    def __next__(self):
        fasta_record = next(self.fasta_reader)
        qual_record = next(self.qual_reader)
        
        name = fasta_record.name
        if name != qual_record.name:
            raise FormatError("The read names in the FASTA and QUAL file do "
                              "not match ({0!r} != {1!r})".format(
                              fasta_record.name, qual_record.name))
        try:
            qualities = convert_qualities(qualread.sequence.split(b' '))
        except KeyError as e:
            raise FormatError("Within read named {0!r}: Found invalid quality "
                              "value {1}".format(fasta_record.name, e))
        
        return self.sequence_class(
            name,
            fasta_record.sequence,
            qualities)
    
    def close(self):
        self.fasta_reader.close()
        self.qual_reader.close()

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
    def __init__(self, path, file_format, **kwargs):
        super(InterleavedFileReader, self).__init__(
            path, 'rb', file_format, **kwargs)
    
    def __next__(self):
        return self.create_record(self.file_format.read_pair(self.fileobj))
    
    def iter_single_end(self, end):
        end -= 1
        for reads in itr(self):
            yield reads[end]
