#!/usr/bin/python

import sys

filename=sys.argv[1] # second word

try:
    f = open(filename)
    seqs={}
    name=''
    for line in f:
        line=line.rstrip()
        if line[0] == '>':
            words=line.split()
            name=words[0][1:] # remove the first character
            seqs[name]=''
        else:
            seqs[name]=seqs[name]+line
    print(seqs)
    f.close()
except IOError:
    print('Error with file - might not exist')