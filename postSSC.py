#!/usr/local/Python/3.6.3/bin/python3

def postSSCparser(inputFile,outputFile):

	import subprocess
	import argparse
	import os, sys, re
	import gzip

	if inputFile.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(inputFile)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(inputFile))
			exit(1)
		else:
			inputFile = inputFile[:-3]

	with open(inputFile, 'r') as sscVars, open(outputFile,'w') as parsedVars:
		parsedVars.write('#Chr\tPos\tType\twt\tmu\tT_DP\tT_AF\tN_DP\tNAF\n')
		for lines in sscVars:
			lines = lines.rstrip()
			if lines[0] != "#":
				linesSep = lines.split()
				Chr = linesSep[0]
				Pos = linesSep[1]
				wt = linesSep[3]
				mu = linesSep[4]
				if len(wt) == 1:
					Type = 'SNV'
				if len(wt) > 1:
					if len(mu) < len(wt):
						Type = 'DEL'
					if len(mu) > len(wt):
						Type = 'INS'
				info = linesSep[7].split(';')
				for bits in info:
					if bits.startswith('T_AF='):
						tAF = bits[5:]
					if bits.startswith('T_DP='):
						tDP = bits[5:]
					if bits.startswith('N_AF='):
						nAF = bits[5:]
					if bits.startswith('T_DP='):
						nDP = bits[5:]
				parsedVars.write(Chr+'\t'+Pos+'\t'+Type+'\t'+wt+'\t'+mu+'\t'+tDP+'\t'+tAF+'\t'+nDP+'\t'+nAF+'\n')

	cmdLine = "gzip {}".format(inputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(inputFile))
		exit(1)
	
	cmdLine = "gzip {}".format(outputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(outputFile))
		exit(1)
