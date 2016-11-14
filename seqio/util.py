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
