# Overview

`seqio` is a python library (3.3+ only) for highly optimized reading/writing of *raw* (*i.e.* unaligned) NGS data in a variety of formats:

* FASTA
* FASTQ
* Interleaved FASTQ
* SAM/BAM/CRAM
* Common compression formats (gz, bz2, xz), which are handled automatically (via [xphyle](https://github.com/jdidion/xphyle))

# Details

`seqio` achieves significant performance gains over similar libraries (e.g. `screed`) in three ways:

1. Implements critical portions of code in C (via Cython).
2. Maintains data as bytes by default, rather than perform expensive string encoding/decoding operations.
3. Uses system-level compression/decompression (via [xphyle](https://github.com/jdidion/xphyle)) when possible.

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

These are all installable via pip:

* Cython 0.25+ (only required at build-time)
* xphyle 1.0+
* pysam (only required for SAM/BAM/CRAM support)

# TODO

* Add FAST5 support (using poretools)
* Add xphyle plugins for sequence-specific compression formats:
    * Treat BAM and CRAM as compression formats using native samtools/sambamba if available
    * fqzcomp, quip, scalce http://www.nature.com.ezproxy.nihlibrary.nih.gov/nmeth/journal/v13/n12/pdf/nmeth.4037.pdf
    * https://github.com/COMBINE-lab/quark
    * https://github.com/rcanovas/libCSAM
* Casava/BCL support might be useful to some; will wait to see if this is requested.
