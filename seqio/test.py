import sys

AT = ord('@')
PLUS = ord('+')

infile = sys.argv[1]
with open(infile, 'rb') as i:
    for name, sequence, name2, qualities in zip(*[i] * 4):
        assert name[0] == AT
        assert name2[0] == PLUS
        name = name[1:-1].decode()
        if len(name2) == 1:
            name2 = None
        else:
            name2 = name2[1:-1].decode()
        sequence = sequence.rstrip()
        qualities = qualities.rstrip()
        
