# -*- coding: utf-8 -*-
"""Utility classes/methods.
"""
class OptionalDependency(object):
    """Subclass property to make classmethod properties possible.
    """
    def __init__(self, name):
        self.name = name
    
    @property
    def lib(self):
        """Loads the python module and replaces itself with that module.
        
        Returns:
            The module
        """
        lib = import_module(self.name)
        self.lib = lib
        lib

class BatchIterator(object):
    def __init__(self, reader, size, max_reads=None):
        self.reader = reader
        self.iterable = enumerate(reader, 1)
        self.size = size
        self.max_reads = max_reads
        self.done = False
        self._empty_batch = [None] * size
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.done:
            raise StopIteration()
        
        try:
            read_index, record = next(self.iterable)
        except StopIteration:
            self.close()
            raise
        
        batch = copy.copy(self._empty_batch)
        batch[0] = record
        batch_index = 1
        max_size = self.size
        if self.max_reads:
            max_size = min(max_size, self.max_reads - read_index + 1)
        
        while batch_index < max_size:
            try:
                read_index, record = next(self.iterable)
                batch[batch_index] = record
                batch_index += 1
            except StopIteration:
                self.close()
                break
        
        if self.max_reads and read_index >= self.max_reads:
            self.close()
        
        if batch_index == self.size:
            return (batch_index, batch)
        else:
            return (batch_index, batch[0:batch_index])
    
    def close(self):
        self.done = True
        self.reader.close()

complement = {
    b'A' : b'T',
    b'C' : b'G',
    b'R' : b'Y',
    b'S' : b'S',
    b'W' : b'W',
    b'K' : b'M',
    b'B' : b'V',
    b'D' : b'H',
    b'N' : b'N'
}
for k,v in list(complement.items()):
    complement[v] = k
    complement[k.lower()] = v.lower()
    complement[v.lower()] = k.lower()

def reverse_complement(seq):
    return b''.join(complement[base] for base in reversed(seq))

def sequence_names_match(r1, r2):
    """
    Check whether the sequences r1 and r2 have identical names, ignoring a
    suffix of '1' or '2'. Some old paired-end reads have names that end in '/1'
    and '/2'. Also, the fastq-dump tool (used for converting SRA files to FASTQ)
    appends a .1 and .2 to paired-end reads if option -I is used.
    """
    name1 = r1.name.split(None, 1)[0]
    name2 = r2.name.split(None, 1)[0]
    if name1[-1:] in '12' and name2[-1:] in '12':
        name1 = name1[:-1]
        name2 = name2[:-1]
    return name1 == name2
