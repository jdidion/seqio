


    @abstractmethod
    def iterate(
            self, record_type: Optional[Type[RecordType]] = DefaultRecord
            ) -> Iterator[RecordType]:
        """Iterate over sequence records.

        Args:
            record_type: Type of record to yield, or DefaultRecord if None.

        Yields:
            Instances of RecordType.
        """
        pass

    def dump(
            self, output_prefix: Optional[str] = None,
            output_type: str = 'file', output_format: str = 'fastq',
            writer_kwargs: Optional[dict] = None, format_kwargs: Optional[dict] = None
            ) -> dict:
        """Stream reads from a remote source to output file(s).

        If `output_type` is 'fifo', 'pv' must be callable or `writer_args` must be a
        dictionary with at least the 'buffer' key and the value being a string
        specifying the program to use for buffering instead of pv.

        Otherwise, output file(s) will be compressed using gzip unless `format_kwargs`
        dict is given with a 'compression' key and value of either False or a different
        compression scheme (e.g. 'gz', 'bz2', or 'xz').

        Args:
            output_prefix: Output file prefix. If None, the accession is used.
            output_type: Type of output ('buffer', 'file', or 'fifo').
            output_format: Format of the output file(s).
            writer_kwargs: Additional keyword arguments to pass to the Writer
                constructor.
            format_kwargs: Additional keyword arguments to pass to the FileFormat
                constructor.

        Returns:
            A dict containing the output file names ('file1' and 'file2'),
            and read_count.
        """
        with self:
            # some readers need to begin iterating before they know if the data is
            # paired
            read_iter = iter(self)
            first_read = next(read_iter)
            read_indexes = (1, 2) if self.paired else (1,)

            file_args = dict(
                (f'file{read_idx}', f'{output_prefix or self.accession}.{read_idx}.fq')
                for read_idx in read_indexes)

            if writer_kwargs is None:
                writer_kwargs = {}
                compression = True
            else:
                compression = writer_kwargs.get('compression', True)
            if compression is True:
                compression = 'gz'
            if compression:
                file_args = dict(
                    (key, f'{name}.{compression}')
                    for key, name in file_args.items())

            writer_kwargs.update(file_args)
            writer_kwargs['compression'] = compression

            writer = get_writer(output_type)(**writer_kwargs)

            with get_file_format(output_format)(writer, **(format_kwargs or {})) as fmt:
                fmt(*first_read)
                for read in read_iter:
                    fmt(*read)

        summary = dict(writer_kwargs)
        summary['accession'] = self.accession
        summary['read_count'] = self.read_count
        return summary


class Writer(metaclass=ABCMeta):
    """Interface for classes that write strings to files.

    Args:
        paired: Whether the reads are paired-end.
    """
    def __init__(self, paired: bool = False):
        self.paired = paired

    @abstractmethod
    def __call__(self, read1_str: str, read2_str: Optional[str] = None):
        """Write strings to a pair of files.

        Args:
            read1_str: The read1 string to write.
            read2_str: The read2 string to write (if paired-end).
        """
        pass

    def close(self) -> None:
        """Close the underlying files.
        """
        pass


def list_writers() -> Tuple[str, ...]:
    return _list_entry_point_names('ngstream.writer')


def get_writer(name: str, dist: Optional[str] = None) -> Callable[..., Protocol]:
    return _get_entry_point(dist, 'ngstream.writer', name)


class Formatter(metaclass=ABCMeta):
    @property
    @abstractmethod
    def paired(self) -> bool:
        pass

    @abstractmethod
    def __call__(self, read1: Record, read2: Optional[Record] = None) -> None:
        """Write a read/pair to the sink.

        Args:
            read1: read1.
            read2: read2 if this is a paired writer, otherwise None.
        """
        pass

    def flush(self) -> None:
        """Flush any buffered output to the sink.
        """
        pass

    def close(self) -> None:
        """Close the sink.
        """
        self.flush()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()


def list_formatters() -> Tuple[str, ...]:
    return _list_entry_point_names('ngstream.formatter')


def get_formatter(name: str, dist: Optional[str] = None) -> Callable[..., Formatter]:
    return _get_entry_point(dist, 'ngstream.formatter', name)



"""Output file formats.
"""
from abc import ABCMeta, abstractmethod
import copy
import os
from typing import List, Optional, cast
from ngstream.api import Formatter, Writer, Record, FastqRecord


class BatchFormatter(Formatter, metaclass=ABCMeta):
    """Wrapper for a string writer (e.g. FifoWriter) that improves performance
    by buffering a set number of reads and sending them as a single call to the
    string writer.

    Args:
        writer: The string writer to wrap. Must be callable with two arguments
            (read1 string, read2 string).
        batch_size: The size of the read buffer.
        lines_per_row: The number of lines used by each read for the specific
            file format (should be passed by the subclass in a
            super().__init__ call).
        linesep: The separator to use between each line (defaults to os.linesep)
    """
    def __init__(
            self, writer: Writer, batch_size: int,
            lines_per_row: int, linesep: str = os.linesep):
        self.writer = writer
        self.batch_size = batch_size
        self.lines_per_row = lines_per_row
        self.bufsize = batch_size * lines_per_row
        self.linesep = linesep
        self._end_records = [self.linesep]
        self.read1_batch = self._create_batch_list()
        if self.paired:
            self.read2_batch = copy.copy(self.read1_batch)
            self._end_records.append(self.linesep)
        self.index = 0

    @property
    def paired(self) -> bool:
        return self.writer.paired

    def _create_batch_list(self) -> List:
        """Create the list to use for buffering reads. Can be overridden, but
        must return a list that is of size ``batch_size * lines_per_row``.
        """
        return [None] * self.bufsize

    def __call__(self, read1: Record, read2: Optional[Record] = None) -> None:
        self.add_to_batch(read1, batch=self.read1_batch, index=self.index)
        if read2:
            self.add_to_batch(read2, batch=self.read2_batch, index=self.index)
        self.index += self.lines_per_row
        if self.index >= self.bufsize:
            self.flush()

    @abstractmethod
    def add_to_batch(self, record: Record, batch: List, index: int) -> None:
        """Add a read to the batch. Must be implemented by a subclass.

        Args:
            record: Read data.
            batch: The list to which data is added.
            index: The item index.
        """
        pass

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def flush(self) -> None:
        """Flush the current read buffers to the underlying string writer.
        """
        def batch_to_str(batch):
            if self.index < self.bufsize:
                return self.linesep.join(batch[0:self.index])
            else:
                return self.linesep.join(batch)
        reads = [batch_to_str(self.read1_batch)]
        if self.paired:
            reads.append(batch_to_str(self.read2_batch))
        self.writer(*reads)
        self.writer(*self._end_records)
        self.index = 0

    def close(self) -> None:
        """Clear the buffers and close the underlying string writer.
        """
        if self.index > 0:
            self.flush()
        self.read1_batch = None
        if self.paired:
            self.read2_batch = None
        self.writer.close()


class FastqFormatter(BatchFormatter):
    """BatchWriter implementation for FASTQ format.
    """
    def __init__(self, writer: Writer, batch_size: int = 1000):
        super().__init__(writer, batch_size, lines_per_row=4)

    def _create_batch_list(self) -> List:
        return [None, None, '+', None] * self.batch_size

    def add_to_batch(
            self, record: Record, batch: List, index: int
            ) -> None:
        batch[index] = '@' + record.name
        batch[index+1] = record.sequence
        if isinstance(record, FastqRecord):
            fastq_record = cast(FastqRecord, record)
            if fastq_record.name2:
                batch[index+2] = f'+{fastq_record.name2}'
        batch[index+3] = record.qualities



"""Writing reads to files.
"""
from io import StringIO
from pathlib import Path
from subprocess import Popen, PIPE
from typing import Tuple, Union, Optional
from xphyle import xopen
from ngstream.api import Writer


class BufferWriter(Writer):
    """StringWriter that writes contents to an in-memory buffer.
    """
    def __init__(self, paired: bool = False):
        super().__init__(paired)
        self._buffers = (StringIO(), StringIO() if paired else None)
        self._values = [None, None]

    @property
    def value(self) -> Tuple[Optional[str], Optional[str]]:
        return (
            self._values[0] or self._buffers[1].getvalue(),
            self._values[1] or (self._buffers[1] if self.paired else None)
        )

    def __call__(self, read1_str: str, read2_str: Optional[str] = None):
        self._buffers[0].write(read1_str)
        if self.paired:
            self._buffers[1].write(read2_str)

    def close(self) -> None:
        self._values[0] = self._buffers[0].getvalue()
        if self.paired:
            self._values[1] = self._buffers[1].getvalue()


class FifoWriter(Writer):
    """StringWriter that opens and writes to a pair of FIFOs in a non-blocking
    way. Each FIFO is written to by opening a subprocess in which stdin is
    piped through pv (with -B option to set max buffer size) to a FIFO.

    Args:
        file1: Path to the read1 FIFO
        file2: Path to the read2 FIFO
        buffer: Command line to program that accepts and buffers between stdin
            and stdout
        kwargs: Additional arguments to pass to Popen
    """
    def __init__(
            self, file1: Path, file2: Optional[Path] = None, buffer='pv -B ',
            compression: Optional[Union[bool, str]] = None, **kwargs):
        super().__init__(file2 is not None)
        self.fifo1 = self._make_fifo(buffer, file1, **kwargs)
        if self.paired:
            self.fifo2 = FifoWriter._make_fifo(buffer, file2, **kwargs)

    @staticmethod
    def _make_fifo(buffer: str, fifo: Path, **kwargs) -> Popen:
        # TODO: use xphyle process instead
        return Popen(
            f'{buffer} > {fifo}', stdin=PIPE, shell=True, universal_newlines=True,
            **kwargs)

    def __call__(self, read1_str: str, read2_str: Optional[str] = None) -> None:
        self.fifo1.stdin.write(read1_str)
        if self.paired:
            self.fifo2.stdin.write(read2_str)

    def close(self) -> None:
        def close_fifo(fifo):
            fifo.stdin.close()
            fifo.terminate()
        close_fifo(self.fifo1)
        if self.paired:
            close_fifo(self.fifo2)


class FileWriter(Writer):
    """StringWriter that opens and writes to a pair of files.

    Args:
        file1: Path to the read1 file.
        file2: Path to the read2 file.
        kwargs: Additional arguments to pass to the ``open`` call.
    """
    def __init__(self, file1: Path, file2: Optional[Path] = None, **kwargs):
        super().__init__(file2 is not None)
        self.file1 = xopen(file1, 'wt', **kwargs)
        if self.paired:
            self.file2 = xopen(file2, 'wt', **kwargs)

    def __call__(self, read1_str: str, read2_str: Optional[str] = None) -> None:
        self.file1.write(read1_str)
        if self.paired:
            self.file2.write(read2_str)

    def close(self) -> None:
        self.file1.close()
        if self.paired:
            self.file2.close()
