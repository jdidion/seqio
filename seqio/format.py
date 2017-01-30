# -*- coding: utf-8 -*-
"""
"""
import textwrap
from xphyle import open_

# some commonly used byte sequences
EMPTY = b''
HASH = b'#'
NEWLINE = b'\n'

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

class TextSequenceFormat(SequenceFormat):
    def __init__(self, sequence_class=Sequence, line_length=None):
        super(TextSequenceFormat, self).__init__(sequence_class)
        self.text_wrapper = None
        if line_length:
            self.text_wrapper = textwrap.TextWrapper(width=line_length)
    
    def _get_sequence(record):
        if self.text_wrapper:
            return self.text_wrapper.fill(
                record.get_full_sequence_str()).encode()
        else:
            return record.full_sequence
    
    def _get_qualities(record):
        if self.delivers_qualities and record.has_qualities:
            return self.text_wrapper.fill(record.get_qualities_str()).encode()
        else:
            return record.qualities
