#!/usr/local/Python/3.6.3/bin/python3


import subprocess
import argparse
import os, sys, re
import operator

#chrom = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']
chrom = ['12']
#'HSV1_GBM_IDT_Probes_B37.bed'
for Chr in chrom:
	posDict = dict()
	with open('tempBed.bed','r') as bedFile:
		for lines in bedFile:
			lines = lines.rstrip()
			linesSplit = lines.split()
			if linesSplit[0] == Chr:
				if linesSplit[1] not in posDict:
					posDict[linesSplit[1]] = linesSplit[2]
				else:
					if int(posDict[linesSplit[1]]) < int(linesSplit[2]):
						posDict[linesSplit[1]] = linesSplit[2]
				posDictSorted = sorted(posDict.items(), key=operator.itemgetter(0))	
	newPosDictSorted = dict()
	for tup in posDictSorted:
		newPosDictSorted[tup[0]] = tup[1]
	print(newPosDictSorted)
	skipper = 0
	for key1 in newPosDictSorted:
		if skipper == 0:
			for key2 in newPosDictSorted:
				if key2 > key1:
					if int(newPosDictSorted[key1]) > int(key2):
						if int(newPosDictSorted[key1]) <= int(newPosDictSorted[key2]):
							print('New setting:', key1, newPosDictSorted[key2])
						elif int(newPosDictSorted[key1]) > int(newPosDictSorted[key2]):
							print('Same setting:', key1, newPosDictSorted[key1])
						skipper = 1
						break
		else:
			print('This is being skipped:', key1)
			skipper = 0
			continue
					#if int(newPosDictSorted[key1]) > int(key2):
					#	print('This is a case!')
