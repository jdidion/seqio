from libcpp.vector cimport vector

cdef class Sequence(object):
    """
    A sequence record has a name and sequence, and optionally
    base qualities and an alternate name. Qualities are encoded
    as ascii(qual+33) by default.
    """
    cdef:
        public str name
        public str name2
        public bytes sequence
        public bytes qualities
        public int length
    
    def __init__(self, str name, bytes sequence, bytes qualities=None, str name2=None):
        self.name = name
        self.name2 = name2
        self._update_sequence(sequence, qualities)
    
    def _update_sequence(self, bytes sequence, bytes qualities=None):
        self.length = len(sequence)
        if qualities and self.length != len(qualities):
            raise FormatError("In read named {0!r}: length of quality sequence "
                "({1}) and length of read ({2}) do not match".format(
                    shorten(name), len(qualities), self.length))
        self.sequence = sequence
        self.qualities = qualities
    
    @property
    def sequence_str(self, **kwargs):
        """Returns the sequence as a string."""
        return self.sequence.decode(**kwargs)
    
    @property
    def has_qualities(self):
        """Whether this sequence has base quality information."""
        return self.qualities is not None:
    
    @property
    def qualities_str(self, **kwargs):
        """Returns qualities as a phred-encoded string."""
        return self.qualities.decode(**kwargs)
    
    def qualities_int(self, int base=33):
        """Returns qualities as a list of integers."""
        assert self.has_qualities
        return list(q - base for q in self.qualities)
    
    def __getitem__(self, key):
        """
        Returns a new Sequence instance with the same name(s) but
        with the sequence and qualities shortened according to the
        given slice.
        """
        return self.__class__(
            self.name,
            self.sequence[key],
            self.qualities[key] if self.has_qualities else None,
            self.name2)

    def __repr__(self):
        qstr = ''
        if self.has_qualities:
            qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
        return '<Sequence(name={0!r}, sequence={1!r}{2})>'.format(
            _shorten(self.name), _shorten(self.sequence), qstr)

    def __len__(self):
        return self.length

    def __richcmp__(self, other, int op):
        """
        Implements == and !=.
        """
        if 2 <= op <= 3:
            eq = (self.name == other.name and
                self.sequence == other.sequence and
                self.qualities == other.qualities)
            if op == 2:
                return eq
            else:
                return not eq
        else:
            raise NotImplementedError()

    def __reduce__(self):
        return (Sequence, (self.name, self.sequence, self.qualities, self.name2))

cdef class ColorspaceSequence(Sequence):
    """
    In colorspace, the first character is the last nucleotide of the primer base
    and the second character encodes the transition from the primer base to the
    first real base of the read.
    """
    cdef public bytes primer
    
    def __init__(self, str name, bytes sequence, bytes qualities=None, bytes primer=None, str name2=None):
        if primer is None:
            assert len(sequence) > 0
            primer = sequence[0:1]
            sequence = sequence[1:]
        
        if not primer in ('A', 'C', 'G', 'T'):
            raise FormatError("Primer base is {0!r} in read {1!r}, but it should be one of A, C, G, T.".format(
                primer, shorten(name)))
        
        super(Colorspace, self).__init__(name, sequence, qualities, name2)
        self.primer = primer

    def __repr__(self):
        qstr = ''
        if self.qualities is not None:
            qstr = ', qualities={0!r}'.format(shorten(self.qualities))
        return '<ColorspaceSequence(name={0!r}, primer={1!r}, sequence={2!r}{3})>'.format(
            shorten(self.name), self.primer, shorten(self.sequence), qstr)

cdef bytes EMPTY = b''

ctypedef struct Edit:
    int start
    int stop
    bytes bases
    bytes quals
    int size_change
    str description

cdef class Mutable(object):
    """
    Sequence mixin that adds edit operations and replaces the __getitem__ method,
    such that the sequence and qualities are updated in-place. Edits are logged
    in the 'edits' list, which enables the original sequence to be reconstructed.
    """
    cdef public vector[Edit] edits
    
    def edit(self, int start=0, int stop=-1, bytes bases=EMPTY,
             bytes qualities=EMPTY, str description=''):
        """
        Modify the current sequence and qualities. The current bases/qualities
        between start (inclusive) and stop (exclusive) are replaced by `bases`,
        and a new Edit operation is generated, added to `edits`, and returned.
        """
        cdef:
            int cur_size = len(self)
            bytes new_sequence = None
            bytes new_qualities = None
            Edit e
        
        if stop < 0:
            stop = cur_size
        elif stop > cur_size:
            raise Exception("Sequence is shorter than edit region {} < {}".format(stop, cur_size))
        
        e = Edit()
        e.start = start
        e.stop = stop
        e.size_change = len(new_bases) - (stop - start)
        e.description = description
        e.bases = self.sequence[start:stop]
        new_sequence = self.sequence[:start] + bases + self.sequence[stop:]
        if self.has_qualities:
            e.quals = self.qualities[start:stop]
            new_qualities = self.qualities[:start] + bases + self.qualities[stop:]
        self._update_sequence(new_sequence, new_qualities)
        self.edits.push_back(edit)
        return edit
        
        # This code might make things slightly faster, at the cost
        # of additional complexity:
        # if start < 0 and stop < 0:
        #     start = 0
        #     stop = cur_size
        #     e.bases = self.sequence
        #     self.sequence = new_bases
        #     if self.has_qualities:
        #         e.quals = self.qualities
        #         self.qualities = new_qualities
        # elif start < 0:
        #     start = 0
        #     e.bases = self.sequence[:stop]
        #     self.sequence = bases + self.sequence[stop:]
        #     if self.has_qualities:
        #         e.quals = self.qualities[:stop]
        #         self.qualities = qualities + self.qualities[stop:]
        # elif stop < 0:
        #     stop = cur_size
        #     e.bases = self.sequence[start:]
        #     self.sequence = self.sequence[:start] + bases
        #     if self.has_qualities:
        #         e.quals = self.qualities[start:]
        #         self.qualities = self.sequence[:start] + qualities
        # else:
    
    def insert(self, int pos, bytes bases, bytes qualities=EMPTY, str description=''):
        """
        Convenience method to insert a sequence directly *before* `pos`.
        """
        return self.edit(pos, pos, bases, qualities, description)
    
    def delete(self, int start=0, int stop=-1, str description=''):
        """
        Convenience method to delete the sequence between `start` (inclusive)
        and `stop` (exclusive).
        """
        return self.edit(start, stop, description=description)
    
    def __getitem__(self, key):
        """
        Generates up to two deletion events - one from the front
        of the read (if `key.start` > 0) and one from the end of
        the read (if `key.stop` < len(self)).
        """
        if key.stop < len(self):
            self.delete(start=key.stop)
        if key.start > 0:
            self.delete(stop=key.start)
        return self

cdef class MutableSequence(Sequence, Mutable):
    """
    Inherits both Sequence and Mutable and adds no additional
    functionality.
    """
    pass

cdef class MutableColorspaceSequence(ColorspaceSequence, Mutable):
    """
    Inherits both ColorspaceSequence and Mutable and adds no additional
    functionality.
    """
    pass

def sra_colorspace_sequence(name, sequence, qualities, name2=None, mutable=False):
    """
    Factory for an SRA colorspace sequence (which has one quality value too many).
    """
    assert qualities is not None and len(qualities) > 0
    if mutable:
        return MutableColorspaceSequence(name, sequence, qualities[1:], name2=name2)
    else:
        return ColorspaceSequence(name, sequence, qualities[1:], name2=name2)
