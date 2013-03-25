#!/usr/bin/env python
# simple script to test the UCT optimizer

from random import gauss, random

p = random()

if p < .5:
	print 'W'
else:
	print 'L'