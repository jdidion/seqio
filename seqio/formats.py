# kate: syntax Python;
# cython: profile=False, emit_code_comments=False
"""Cython implementations of sequence formats.
"""
from xphyle import open_

class ClassProperty(property):
    """Subclass property to make classmethod properties possible.
    """
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()

class SequenceFormat(object):
    def open(self, path, mode, **kwargs):
        return open_(path, mode, **kwargs)

class SAM(SequenceFormat):
    """SAM/BAM/CRAM format files. Paired-end files must be name-sorted. Does not
    support secondary/supplementary reads.
    """
    name = 'sam'
    aliases = ('sam', 'bam', 'cram')
    readable = True
    writeable = True
    delivers_qualities = True
    _lib = None
    
    @ClassProperty
    @classmethod
    def lib(self):
        """Caches and returns the python module assocated with this file format.
        
        Returns:
            The module
        
        Raises:
            CompressionError if the module cannot be imported.
        """
        if not SAM._lib:
            SAM._lib = import_module('pysam')
        return SAM._lib
    
    def open(path, mode, **kwargs):
        return self.lib.AlignmentFile(path, mode, **kwargs)
    
    
    
    def __init__(self, path, sequence_class=Sequence, **kwargs):
        super(SAMReader, self).__init__(path, **kwargs)
        self.sequence_class = sequence_class
    
    def __iter__(self):
        
    
    def _as_sequence(self, read):
        return self.sequence_class(
            read.query_name,
            read.query_sequence,
            ''.join(chr(33 + q) for q in read.query_qualities))

class SingleEndSAMReader(SAMReader):
    def _iter(self, sam):
        for read in sam:
            yield self._as_sequence(read)

class Read1SingleEndSAMReader(SAMReader):
    def _iter(self, sam):
        for read in sam:
            if read.is_read1:
                yield self._as_sequence(read)

class Read2SingleEndSAMReader(SAMReader):
    def _iter(self, sam):
        for read in sam:
            if read.is_read2:
                yield self._as_sequence(read)

class PairedEndSAMReader(SAMReader):
    def _iter(self, sam):
        for reads in zip(sam, sam):
            if reads[0].query_name != reads[1].query_name:
                raise Exception(
                    "Consecutive reads {}, {} in paired-end SAM/BAM file do not "
                    "have the same name; make sure your file is name-sorted and "
                    "does not contain any secondary/supplementary alignments.",
                    reads[0].query_name, reads[1].query_name)
            
            if reads[0].is_read1:
                assert reads[1].is_read2
            else:
                assert reads[1].is_read1
                reads = (reads[1], reads[0])
            
            yield tuple(self._as_sequence(r) for r in reads)
