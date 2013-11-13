#!/usr/bin/env python

# From: http://www.michael-noll.com/tutorials/writing-an-hadoop-mapreduce-program-in-python/

import sys


#This function truncates a number f to n decimals
def trunc(f, n):
    return str(f)[:len('%.*f' % (n, f))]

# input comes from STDIN (standard input)
for line in sys.stdin:
    #sys.stdin
    # remove leading and trailing whitespace
    line = line.strip()
    # split the line into words
    line = line.strip() # Need to remove trailing '\n'
    entries = line.split('\t')
    out_str = str(trunc(float(entries[0]),1)) + str(',') + str(float(trunc(float(entries[0]),1))+0.1) + str(',') + str(trunc(float(entries[1]),1)) + str(',') + str(float(trunc(float(entries[1]),1))+0.1) + str('\t') + str(1)

    print out_str

