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

def complement(seq):
    return b''.join(complement[base] for base in seq)

def reverse_complement(seq):
    return b''.join(complement[base] for base in reversed(seq))

# http://nbviewer.jupyter.org/urls/gist.github.com/Jorge-C/d51b48d3e18897c46ea2/raw/73d7e11e4b72d6ba90e0021931afa230e63031e9/cython+sequences.ipynb?create=1
cdef class CySeq3(object):
    """Translation using a LUT, cythonized, uses numpy but str as input"""
    arr = np.array([  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                     84,  86,  71,  72,   0,   0,  67,  68,   0,   0,  77,   0,  75,
                     78,   0,   0,   0,  89,  83,  65,   0,  66,  87,   0,  82,   0,
                      0,   0,   0,   0,   0,   0, 116, 118, 103, 104,   0,   0,  99,
                    100,   0,   0, 109,   0, 107, 110,   0,   0,   0, 121, 115,  97,
                      0,  98, 119,   0, 114,   0,   0,   0,   0,   0,   0], dtype=np.uint8)
    cdef public cnp.ndarray seq
    cdef int length
    
    def __cinit__(self, str seq):
        self.seq = np.array(seq, dtype='c').view(np.uint8)
        self.length = len(seq)


    @cython.boundscheck(False)
    cdef _rc(self):
        cdef cnp.ndarray[cnp.uint8_t, ndim=1] arr = self.arr
        cdef cnp.ndarray[cnp.uint8_t, ndim=1] seq = self.seq
        cdef cnp.ndarray[cnp.uint8_t, ndim=1] result = np.empty_like(seq)
        cdef int i
        for i in reversed(range(self.length)):
            result[-i-1] = arr[seq[i]]
        return result.view('c')

    def rc(self):
        return self._rc()

def qual2prob(qual):
    """Convert a phred-scale integer to a probability.
    """
    return 10 ** (-q / 10)

def prob2qual(prob):
    """Convert a probability to a phred-scale integer.
    """
    if p < 0:
        raise ValueError("Probability must be >= 0")
    return round(-10 * math.log10(prob))

QualityValue = namedtuple('QualityValue', ('int', 'prob', 'asc'))

def guess_quality_type(q):
    # val as an integer phred score
    if isinstance(q, int):
        return 'int'
    # val is a probability
    if isinstance(q, float):
        return 'prob'
    # val is an ascii-encoded phred score
    if isinstance(q, str):
        return 'asc'
    if isinstance(q, bytes):
        return 'asc'
    raise ValueError("Cannot guess quailty type: {}".format(q))

class QualityConversion(object):
    """Converts between phred-scaled int scores, ascii-encoded phred scores,
    and probabilities.
    
    Args:
        base: The base for ascii-to-int conversion. Typically 33 unless
            older Illumina qualities are used, in which case the base
            should be 64.
        max_asc: The max ascii value to allow. Must be <= 255.
    """
    def __init__(self, base=33, max_asc=255):
        if max_asc is None or max_asc < 0 or max_asc > 255:
            raise ValueError("'max_asc' must be 0 < i <= 255")
        #if base not in (33, 64):
            #this should be a warning
            #pass
        self.base = base
        self.max_i = max_asc - base
        self.int_table = {}
        self.asc_table = {}

    def convert_iter(quals, dest, src=None):
        """Convert qualities from one format to another.
        
        Args:
            quals: An iterable of qualities
            dest: The destination format ('int', 'asc', 'prob')
            src: The source format -- one of the above or None if it should be
                guessed
            base: The base for ascii-to-int conversion
        
        Returns:
            An iterable of converted qualities
        """
        if not quals:
            return ()
        if src is None:
            src = guess_quality_type(quals[0])
        if src == prob:
            quals = (prob2qual(q) for q in quals)
            src = 'int'
        return (self.convert(q, dest, src) for q in quals)
    
    def convert(self, qual, dest, src=None):
        """Convert a quality score from one format to another.
        
        Args:
            qual: A quality score
            dest: The destination format ('int', 'asc', 'prob')
            src: The source format -- one of the above or None if it should be
                guessed
        
        Returns:
            The converted quality
        """
        if src is None:
            src = guess_quality_type(qual)
        
        a = None
        
        if src == 'prob':
            p = qual
            if p == 0:
                i = self.max_i
            else:
                i = prob2qual(p)
            src = 'int'
        
        if src == 'int':
            if qual in self.int_table:
                return getattr(self.int_table[qual], dest)
            i = qual
        else:
            if qual in self.asc_table:
                return getattr(self.asc_table[qual], dest)
            i = ord(qual) - self.base
            a = qual
        
        if i < 0:
            raise ValueError("Phred score must be >= 0")
        
        i = min(i, self.max_i)
        p = qual2prob(i)
        if a is None:
            a = chr(i + self.base).encode()
        
        record = QualityValue(i, p, a)
        self.int_table[i] = record
        self.asc_table[a] = record
        
        return getattr(record, dest)
