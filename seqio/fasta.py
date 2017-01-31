ARROW = b'>'

class Fasta(TextSequenceFormat):
    name = 'fasta'
    aliases = ('fa',)
    delivers_qualities = False
    
    def __init__(self, linesep=EMPTY, **kwargs):
        super(Fasta, self).__init__(**kwargs)
        self.linesep = linesep
    
    def read_record(self, fileobj):
        while True:
            header = next(fileobj).rstrip()
            if header[0] != ARROW:
                raise FormatError("Expected '>' at beginning of FASTA record")
            elif header[0] == HASH:
                continue
            break
        
        seq = []
        while True:
            if fileobj.peek(1) in (ARROW, EMPTY):
                break
            line = next(fileobj).rstrip()
            if line:
                seq.append(line)
        
        return self.sequence_class(
            name=header[1:], sequence=self.linesep.join(seq))
    
    def format_record(self, record):
        if self.text_wrapper:
            sequence = self.text_wrapper.fill(
                record.get_sequence_str()).encode()
        else:
            sequence = record.sequence
        return EMPTY.join((ARROW, record.name, NEWLINE, sequence, NEWLINE))

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
        self.qual_reader = SingleFileReader(qualfile, Fasta(linesep=b' '), **kwargs)
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

def open(self, *files: FileListArg, mode: str = 'rb', qualities: bool = None,
         format_args: dict = None, io_args: dict = None) -> FileSeqIO:
    if len(files) > 1:
        if qualities is False:
            raise ValueError("More than one file given for FASTA")
        qualities = True
    elif qualities:
        raise ValueError("Two files required for FASTQUAL")
    if qualities:
        klass = FastaQualReader if 'r' in mode else FastaQualWriter
        return klass(*files, mode=mode, format_args=format_args, io_args=io_args)
    else:
        klass = FastaQualReader if 'r' in mode else FastaQualWriter
        file_format = Fastq(**(format_args or {}))
        return klass(files[0], file_format, mode=mode, **(io_args or {}))
