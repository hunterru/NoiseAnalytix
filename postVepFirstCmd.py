#!/usr/local/Python/3.6.3/bin/python3

def postVepParserFirst(inputFile,outputFile,workingDir):

	import subprocess
	import argparse
	import os, sys, re

	tempOutputFile = workingDir+'/tempFirstParser.tsv'
	if inputFile.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(inputFile)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(inputFile))
			exit(1)
		else:
			inputFile = inputFile[:-3]

	with open(inputFile,'r') as vepVars1, open(outputFile,'w') as vepOut:
		pos0sep = '0'
		for line1 in vepVars1:
			line1 = line1.rstrip()
			if line1[0] != "#":
				line1Sep = line1.split()
				if pos0sep == line1Sep[0]:			
					continue
				else:
					cmdLine = "grep '"+'{}'.format(line1Sep[0])+"' "+inputFile+" > "+tempOutputFile
					cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					if cmdLine_run.returncode != 0:
						print("FAILED! Unable to use grep command in postVepFirstCmd.py")
						exit(1)
				pos0sep = line1Sep[0]			

				with open(tempOutputFile,'r') as select:
					sEnst = 0
					for lines in select:
						lines = lines.rstrip()
						linesSep = lines.split()
						Enst = linesSep[4]
						found1 = re.search(r"ENST0+",Enst)
						cEnst = int(Enst[len(found1.group(0)):])
						if cEnst > sEnst:
							sEnst = cEnst
							out = linesSep
					impactSplit = out[6].split(',')
					for s in impactSplit:
						if (s == 'stop_gained') or (s == 'missense_variant'):
							vepOut.write(lines+'\n')
							break

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

	cmdLine = "rm {}".format(tempOutputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to remove {}".format(tempOutFile))
		exit(1)
