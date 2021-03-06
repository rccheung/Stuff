#!/usr/bin/env python

from operator import itemgetter
import sys

current_key = None
current_count = 1

# input comes from STDIN
for line in sys.stdin:
    # remove leading and trailing whitespace
    line = line.strip()

    # parse the input we got from mapper.py
    key, count = line.split('\t'); count = int(count)
    try:
        count = int(count)
    except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
        continue

    # this IF-switch only works because Hadoop sorts map output
    # by key (here: word) before it is passed to the reducer
    if current_key == key:
        current_count += count
    else:       #new xlo
        if current_key:
            keys = current_key.split(',')
            xlo = float(keys[0]); xhi = float(keys[1]); ylo = float(keys[2]); yhi = float(keys[3]); count = int(count)
            print '%.1f,%.1f,%.1f,%.1f,%d' % (xlo,xhi,ylo,yhi,current_count)
        current_count = count
        current_key = key  #reset count

# do not forget to output the last word if needed!
if current_key == key:
    keys = current_key.split(',')
    xlo = float(keys[0]); xhi = float(keys[1]); ylo = float(keys[2]); yhi = float(keys[3]); count = int(count)            
    print '%.1f,%.1f,%.1f,%.1f,%d' %(xlo,xhi,ylo,yhi,current_count)


