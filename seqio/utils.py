# -*- coding: utf-8 -*-
"""Utility classes/methods.
"""

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
