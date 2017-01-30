#!/usr/bin/env python
from xphyle import open_
from random import choices

def gen_seq(length, qualities=True):
    seq = ''.join(choices('ACGT', k=length))
    if qualities:
        return (seq, 'I' * length)
    else:
        return seq

def gen_fastq(n, start, id_template, length):
    return [
        (id_template.format(i),) + gen_seq(length)
        for i in range(start, start+n)
    ]

def gen_fastq_series(m, name_template, n, id_template, length):
    files = {}
    for pair in (1, 2):
        i = 1
        for idx in range(65, 65+m):
            name = name_template.format(chr(idx), pair)
            files[name] = gen_fastq(n, i, id_template, length)
            i += n
    return files

fastq = gen_fastq_series(2, 'test{}.{}.fq.gz', 2, 'rec{}', 8)

if __name__ == '__main__':
    for fname, content in fastq.items():
        with open_(fname, 'w') as o:
            o.write('\n'.join("{}\n{}\n+\n{}".format(*rec) for rec in content))
