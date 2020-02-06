#!/usr/bin/env python

from BioInfoUtils import *

fn = raw_input("file? ")
fp = open(fn, 'r')
text = fp.readline().rstrip()

print translate(text)
