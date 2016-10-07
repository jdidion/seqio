# Overview

`seqio` is a python library (3.3+ only) for highly optimized reading/writing of common NGS file formats. Specifically, `seqio` supports:

* FASTA (r+w)
* FASTQ (r+w)
* Interleaved FASTQ (r+w)
* SAM/BAM/CRAM (read-only)
* Common compression formats (gz, bz, lz), which are handled automatically (via `gemisch`)

# Details

`seqio` achieves significant performance gains over similar libraries (e.g. `screed`) in three ways:

1. Implements critical portions of code in C (via cython).
2. Maintains sequences and qualities as byte strings by default, rather than perform expensive encoding/decoding operations.
3. Uses system-level compression/decompression (via `gemisch`) when possible.

The generated sequences can be either immutable or mutable. With immutable sequences, slicing always returns a new sequence, and modifications to the sequence and qualities are not supported. With mutable sequences, modifications are applied in-place and each modification is recorded as an `Edit`.

# Usage

```python
# Here we will filter out reads with low mean quality
# or which contain an adapter sequence, and write the
# retained reads to an interleaved, compressed fastq file
from seqio import reader, writer

fq = ('reads{}.fq.gz'.format(i) for i in (1,2))
adapter = b'AGATCGAAT'
writer = writer('output.fq.gz', interleaved=True)

for read1, read2 in reader(*fq):
    # qualities can also be retrieved as a list of
    # phred-scale integers
    mean_qual = sum(read1.qualities_phred(base=33)) / len(read1)
    # the following is equivalent and a bit faster:
    from seqio.util import mean_quality
    mean_qual = mean_quality(read1)
    
    # for most operations (e.g. comparisons), byte array
    # sequence and qualities can be used directly
    has_adapter = read1.sequence.startswith(adapter)
    # the following is equivalent and will take care of
    # string conversion for us, but is slightly slower
    from seqio.util import has_prefix
    has_adapter = has_prefix(read1, 'AGATCGAAT')
    
    # names are stored as unicode strings
    print(read1.name)
    # sequences and qualities are stored as byte arrays,
    # so explicit conversion to string is required
    print("{}\n{}\n\n{}\n{}".format(
        read1.sequence_str, read1.qualities_str,
        read2.sequence_str, read2.qualities_str))
    print("Mean read1 quality: {:.3d}".format(mean_qual))
    print("Contains adapter? {}".format(has_adapter))
    
    if mean_qual >= 30 and not has_adapter:
        writer.write(read1, read2)
```

# Dependencies

* Cython 0.25+
* gemisch
* pysam (required for SAM/BAM/CRAM support)
