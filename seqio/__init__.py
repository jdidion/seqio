# -*- coding: utf-8 -*-
"""
"""
from importlib import import_module
from seqio.types import PathOrFile, ModeArg
from xphyle.utils import is_iterable

class Formats(object):
    """Container for file format instances.
    """
    def __init__(self):
        self.formats = {}
    
    def register(name: str, mod) -> bool:
        """Register a file format.
        
        Args:
            name: The format name.
            mod: The format implementation (a module or object).
        
        Returns:
            True if the format already exists and was replaced, otherwise False.
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

FilesArg = Union[PathSpec, PathOrFile, Iterable[PathOrFile]]
"""A path, iterable of paths, or :class:`xphyle.paths.PathSpec`."""

def open(
        files1: FilesArg, files2: FilesArg = None, mode: ModeArg = 'rb',
        file_format: FormatArg = None, **kwargs) -> SequenceReader:
    """Open a sequence file reader/writer.
    
    Args:
        files1, files2: The file(s) to open. See `Notes`.
        mode: The mode for reading/writing data. Defaults to 'rb'.
        file_format: The file format, or `None` to auto-detect. All files must
            be of the same format.
        kwargs: Format-specific arguments.
    
    Notes:
        `files` can be any of the following:
        
        * A single file (path string or file-like object)
        * A glob (e.g. foo*.1.fq.gz)
        * A tuple of files (to read successively from multiple files)
        * A :class:`xphyle.paths.PathSpec`, which describes how to locate files;
        must satisfy the requirements for the specific file format.
        
        `mode` must be a valid :class:`xphyle.types.FileMode`. When the mode is
        binary, data are read as bytes and may be optionally converted to
        strings by calling the appropriate methods (e.g. name_str,
        sequence_str).
    
    Returns:
        A :class:`SequenceReader`.
    """
    if isinstance(mode, str):
        mode = FileMode(mode)
    
    if not file_format:
        if is_iterable(files1):
            files1 = tuple(files1)
            if files2:
                files2 = tuple(files2)
            test_item = files1[0]
        else:
            test_item = files1
        
        if isinstance(test_item, PathSpec):
            file_format = FORMATS.guess_from_path_spec(test_item, mode.readable)
        else:
            file_format = FORMATS.guess_from_file(test_item, mode.readable)
        
        if file_format is None:
            raise ValueError("Cannot guess file format")
    
    elif isinstance(file_format, str):
        # get file_format from format name
        file_format = get_format(file_format)
    
    files = (files1, files2) if files2 else (files1,)
    return file_format.open(*files, mode=mode, **kwargs)
