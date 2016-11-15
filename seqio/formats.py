# -*- coding: utf-8 -*-
"""
"""
import textwrap
from xphyle import open_

class SequenceFormat(object):
    def __init__(self, sequence_class=Sequence):
        self.sequence_class = sequence_class

    def open(self, path, mode, **kwargs):
        return open_(path, mode, **kwargs)

    def read_record(self, fileobj):
        raise NotImplemented()
    
    def read_pair(self, fileobj):
        return (self.read_record(fileobj), self.read_record(fileobj))
    
    def _create_record(self, *args, **kwargs):
        return self.sequence_class(*args, **kwargs)

    def format_record(self, record):
        raise NotImplemented()
    
    def format_pair(self, read1, read2):
        return (self.format_record(read1), self.format_record(read2))

class FastA(SequenceFormat):
    name = 'fasta'
    aliases = ('fa',)
    delivers_qualities = False
    
    def __init__(self, line_length=None):
        self.text_wrapper = None
        if line_length:
            self.text_wrapper = textwrap.TextWrapper(width=line_length)
    
    def read_record(self, fileobj):
        while True:
            header = next(fileobj).rstrip()
            if header[0] != b'>':
                raise FormatError("Expected '>' at beginning of FASTA record")
            elif header[0] == '#':
                continue
            break
        
        seq = []
        while True:
            if fileobj.peek(1) in (b'>', b''):
                break
            line = next(fileobj).rstrip()
            if line:
                seq.append(line)
        
        return self.sequence_class(name=header[1:], sequence=b''.join(seq))
    
    def format_record(self, record):
        if self.text_wrapper:
            sequence = self.text_wrapper.fill(record.sequence_str).encode()
        else:
            sequence = record.sequence
        return b''.join((b'>', record.name, b'\n', sequence, b'\n'))

class FastQ(SequenceFormat):
    name = 'fastq'
    aliases = ('fq',)
    delivers_qualities = True
    
    def read_record(self, fileobj):
        lines = [next(fileobj).rstrip() for i in range(4)]
        if lines[0][0] != b'@' or lines[2][0] != b'+':
            raise FormatError(
                "FASTQ record is formatted incorrectly: {}".format(
                b''.join(lines)))
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
            name=name, name2=name2,
            sequence=lines[1],
            qualities=lines[3])
    
    def format_record(self, record):
        return b''.join((
            b'@', record.name, b'\n',
            record.sequence, b'\n+',
            record.name2, b'\n',
            record.qualities, b'\n'
        ))

class SAM(SequenceFormat):
    """SAM/BAM/CRAM format files. Paired-end files must be name-sorted. Does not
    support secondary/supplementary reads.
    """
    name = 'sam'
    aliases = ('bam', 'cram')
    delivers_qualities = True
    lib = OptionalDependency('pysam')
    
    def open(path, mode, **kwargs):
        return self.lib.AlignmentFile(path, mode, **kwargs)
    
    def read_record(self, fileobj):
        for record in fileobj:
            if record.is_secondary or record.is_supplementary:
                continue
            return self._create_record(record)
    
    def read_pair(self, fileobj):
        read1 = read2 = None
        for record in fileobj:
            if record.is_secondary or record.is_supplementary:
                continue
            elif record.is_read1 and not read1:
                read1 = self._create_record(record)
            elif record.is_read2 and not read2:
                read2 = self._create_record(record)
            else:
                raise FormatError(
                    "Expected exactly one read1 and one read2 for read pair "
                    "{}".format((read1 or read2).name))
            if read1 and read2:
                break
        if read1.name != read2.name:
            raise FormatError(
                "Consecutive reads {}, {} in paired-end SAM/BAM file do "
                "not have the same name; make sure your file is name-sorted.",
                reads[0].name, reads[1].name)
        return (read1, read2)
    
    def _create_record(self, record):
        return self.sequence_class(
            name=record.query_name,
            sequence=record.query_sequence,
            qualities=''.join(chr(33 + q) for q in record.query_qualities))
    
    def format_record(self, record):
        record = pysam.AlignedSegment()
        record.query_name = record.name
        record.query_sequence = record.sequence_str
        record.flag = 4
        record.query_qualities = pysam.qualitystring_to_array(record.quality_str)
        return record
    
    def format_pair(self, read1, read2):
        record1 = format_record(read1)
        record1.flag = 77
        record2 = format_record(read2)
        record2.flag = 141
        return (record1, record2)