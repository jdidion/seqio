# -*- coding: utf-8 -*-
"""
"""
from importlib import import_module
from seqio.types import FileListArg, PathOrFile
from xphyle.utils import is_iterable

class Formats(object):
    """Container for file format instances.
    """
    def __init__(self):
        self.formats = {}
    
    def register(name: str, mod):
        """Register a file format.
        
        Args:
            name: The format name
            mod: The format implementation (a module or object)
        
        Returns:
            True if the format already exists and was replaced, otherwise False
        """
        exists = name in self.formats
        self.formats[name] = mod
        return exists
    
    def guess_from_file(path_or_file: PathOrFile):
        pass

    def guess_from_path_spec(path_spec):
        pass

FORMATS = Formats()

# register known formats
for fmt in ('fasta', 'fastq', 'sam'):
    mod = import_module('seqio.{}'.format(fmt))
    FORMATS.register(fmt, mod)

def reader(*files: FileListArg, file_format: Union[str, FileFormat] = None,
           mode: str = 'b', **kwargs) -> Reader:
    """Open a sequence file reader.
    
    Args:
        files: The file(s) to open.
        mode: The mode for reading data; 'b' = data are read as bytes and
            may be optionally converted to strings by calling the appropriate
            methods (e.g. name_str, sequence_str); 't' = read data as text.
        file_format: The file format, or None to auto-detect. All files must be
            of the same format.
        kwargs: File type-specific parameters.
    
    'files' can be any of the following:
    * A single file (path string or file-like object)
    * Two files (e.g. for paired-end fastq)
    * A tuple of files (to read successively from multiple files)
    * Two tuples of files (to read paired-end data successively from multiple
      files)
    * An xphyle.paths.PathSpec, which describes how to locate files; must
      satisfy the requirements for the specific file format.
    """
    # get file_format from format name
    if file_format and isinstance(file_format, str):
        file_format = get_format(file_format)
    
    # otherwise detect from file name
    if file_format is None:
        if is_iterable(files[0]):
            files = [tuple(f) for f in files]
            test_item = files[0][0]
        else:
            test_item = files[0]
        
        if isinstance(test_item, PathSpec):
            file_format = FORMATS.guess_from_path_spec(test_item)
        else:
            file_format = FORMATS.guess_from_file(test_item)
    
    if file_format is None:
        raise ValueError("Cannot guess file format")
    
    return file_format.open_reader(*files, mode=mode, **kwargs)
