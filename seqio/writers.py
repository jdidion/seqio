from seqio.base import *

class SequenceWriter(FileSeqIO):
    paired = False
    
    def __init__(self, path, file_format, **kwargs):
        super(SequenceWriter, self).__init__(
            path, 'wb', file_format, **kwargs)
    
    def write(self, record):
        self.fileobj.write(self.file_format.format_record(record))

class PairedFileWriter(SeqIO):
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
