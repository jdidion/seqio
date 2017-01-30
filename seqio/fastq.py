from seqio.format import TextSequenceFormat
from seqio.io import (
    SingleFileReader, PairedFileReader, InterleavedFileReader,
    SingleFileWriter, PairedFileWriter, InterleavedFileWriter)
from seqio.types import FileListArg

AT = b'@'
PLUS = b'+'

class Fastq(TextSequenceFormat):
    name = 'fastq'
    aliases = ('fq',)
    delivers_qualities = True
    
    def __init__(self, write_name2: bool = False, **kwargs):
        super(Fastq, self).__init__(self, **kwargs)
        self.write_name2 = write_name2
    
    def read_record(self, fileobj):
        lines = [next(fileobj).rstrip() for i in range(4)]
        if lines[0][0] != AT or lines[2][0] != PLUS:
            raise FormatError(
                "FASTQ record is formatted incorrectly: {}".format(
                EMPTY.join(lines)))
        name = lines[0][1:]
        name2 = lines[2][1:]
        if name2 and name != name2:
            raise FormatError(
                "Sequence descriptions in the FASTQ file don't match "
                "({0!r} != {1!r}).\n"
                "The second sequence description must be either empty "
                "or equal to the first description.".format(
                name, name2))
        return self._create_record(
            name=name,
            sequence=lines[1],
            qualities=lines[3])
    
    def format_record(self, record):
        return EMPTY.join((
            AT, record.name, NEWLINE,
            record.sequence, NEWLINE,
            PLUS, record.name if self.write_name2 else EMPTY, NEWLINE,
            record.qualities, NEWLINE
        ))

def open_reader(self, *files: FileListArg, mode: str = 'b', interleaved=False,
                paired=None, format_args: dict = None, io_args: dict = None
               ) -> Reader:
    file_format = Fastq(**(format_args or {}))
    if len(files) == 1:
        if interleaved:
            reader_class = InterleavedFileReader
        elif paired:
            raise ValueError("Two files required for paired-end FASTQ")
        else:
            reader_class = SingleEndFileReader
    elif paired is False:
        raise ValueError("Two files required for paired-end FASTQ")
    elif interleaved:
        raise ValueError("Interleaved FASTQ must be a single file")
    else:
        reader_class = PairedEndFileReader
    return reader_class(*files, file_format, mode=mode, **(io_args or {}))
