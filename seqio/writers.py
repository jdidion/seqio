from seqio.base import *

class SingleWriter(object):
    paired = False

class SequenceWriter(FileSeqIO):
    def __init__(self, path, file_format, **kwargs):
        super(SequenceWriter, self).__init__(
            path, 'wb', file_format, **kwargs)
    
    def write(self, record):
        self.fileobj.write(self.file_format.format_record(record))

class FastaQualWriter(SeqIO, SingleWriter):
    """Writer for reads that are stored in .(CS)FASTA and .QUAL files.
    
    Args:
        fastafile: path or file-like object of sequences in FASTA format
        qualfile: path or file-like object of qualities in FASTA format
        kwargs: Additional arguments passed to the file open method
    """
    delivers_qualities = True
    
    def __init__(self, fastafile, qualfile, **kwargs):
        fasta = Fasta()
        self.fasta_reader = SingleFileReader(fastafile, fasta, **kwargs)
        self.qual_reader = FastaReader(qualfile, fasta, **kwargs)
    
    def __next__(self):
        fasta_record = next(self.fasta_reader)
        qual_record = next(self.qual_reader)
        
        name = fasta_record.name
        if name != qual_record.name:
            raise FormatError("The read names in the FASTA and QUAL file do "
                              "not match ({0!r} != {1!r})".format(
                              fasta_record.name, qual_record.name))
        try:
            qualities = Qualities(list(
                int(q) for q in qualread.sequence.split(b' ')), 'int')
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
