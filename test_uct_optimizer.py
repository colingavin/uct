#!/usr/bin/env python
# simple script to test the UCT optimizer

from random import gauss, random
import sys

p = random()

k = float(sys.argv[2])

if k <= 500:
    prob = -k*(k - 500.0)/(62500.0)
else:
    prob = -(k - 500.0)*(k - 1000.0)/(2*62500.0)

#print k, prob

if p < prob:
	print 'W'
else:
	print 'L'
