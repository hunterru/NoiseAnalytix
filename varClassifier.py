#!/usr/local/Python/3.6.3/bin/python3

def postVepVafClass(inputFile,outputFile):

	import subprocess
	import argparse
	import os, sys, re

	with open(inputFile,'r') as VARS, open(outputFile,'w') as outVARS:
		outVARS.write('Chr\tPos\tWT\tMUT\tGENE\tcDesc\tpDesc\tDP\tvarCount\tpVAF\tIMPACT\tSIFT\tPoly\tClinVar\n')
		for lines in VARS:
			lines = lines.rstrip()
			linesSplit = lines.split('\t')
			if linesSplit[0] != 'Chr':
				if float(linesSplit[11]) < 0.05 and float(linesSplit[12]) > 0.80:
					outVARS.write(lines+'\n')
