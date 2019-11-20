#!/usr/local/Python/3.6.3/bin/python3

def commonVars(inputFile1,inputFile2,outputFile):

	import subprocess
	import argparse
	import os, sys, re
	import gzip

	if inputFile1.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(inputFile1)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(inputFile))
			exit(1)
		else:
			inputFile1 = inputFile1[:-3]

	if inputFile2.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(inputFile2)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(inputFile))
			exit(1)
		else:
			inputFile2 = inputFile2[:-3]

	with open(inputFile1,'r') as VARS1, open(outputFile,'w') as commonOut:
		commonOut.write('#Chr\tPos\tType\twt\tmu\tT_DP\tT_AF\tN_DP\tNAF\n')
		counter = 0
		for lines1 in VARS1:
			lines1 = lines1.rstrip()
			if lines1[0] != "#":
				lines1sep = lines1.split()
				Chr = lines1sep[0]
				Pos = lines1sep[1]
				Type = lines1sep[2]
				mu = lines1sep[4]
				with open(inputFile2,'r') as VARS2:
					for lines2 in VARS2:
						lines2 = lines2.rstrip()
						if lines2[0] != "#":
							lines2sep = lines2.split()
							if lines2sep[0] == Chr:
								if lines2sep[1] == Pos:
									if lines2sep[2] == Type:
										if lines2sep[4] == mu:
											commonOut.write(lines1+'\n')
											counter += 1
											break
	
	print('Number of common variants: {}'.format(counter))
	
	cmdLine = "gzip {}".format(outputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(outputFile))
		exit(1)
