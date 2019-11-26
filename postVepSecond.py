#!/usr/local/Python/3.6.3/bin/python3

def postVepParserSecond(inputFile,outputFile):

	import subprocess
	import argparse
	import os, sys, re

	if inputFile.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(inputFile)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(inputFile))
			exit(1)
		else:
			inputFile = inputFile[:-3]

	with open(inputFile,'r') as vepVars, open(outputFile,'w') as dataOut:
		dataOut.write('Chr\tPos\tWT\tMUT\tENSG\tENST\tGENE\tTYPE\tcDesc\tpDesc\tIMPACT\tSIFT\tPoly\tClinVar\n')
		for line in vepVars:
			line = line.rstrip()
			if line != '#':
				lineSep = line.split()
				pos0sep = lineSep[0].split('_')
				Chr = pos0sep[0]
				Pos = pos0sep[1]
				var = pos0sep[2].split('/')
				wt = var[0]
				mu = var[1]
				ENSG = lineSep[3]
				ENST = lineSep[4]
				cDesc = 'c.'+lineSep[7]+wt+'>'+mu
				aa = lineSep[10].split('/')
				pDesc = 'p.'+aa[0]+lineSep[9]+aa[1]
				exParse = lineSep[13].split(';')
				clinVarOut = ''
				for extras in exParse:
					if extras.startswith('SYMBOL='):
						geneName = extras[7:]
					if extras.startswith('SIFT='):
						temp = extras.split('(')
						sift = temp[1][:-1]
					if extras.startswith('PolyPhen='):						
						temp = extras.split('(')
						poly = temp[1][:-1]
					if extras.startswith('CLIN_SIG='):
						temp = extras.split('=')
						clinVar = temp[1].split(',')
						for desc in clinVar:
							clinVarOut = desc
					if extras.startswith('IMPACT='):
						impact = extras[7:]
				typeSplit = lineSep[6].split(',')
				for types in typeSplit:
					if typeSplit[0] == 'stop_gained':
						outData = Chr+'\t'+Pos+'\t'+wt+'\t'+mu+'\t'+ENST+'\t'+ENSG+'\t'+geneName+'\t'+types+'\t'+cDesc+'\t'+pDesc+'\t'+impact
						dataOut.write(outData+'\n')
						break
					elif types == 'missense_variant':
						if (float(sift) <= 0.05) and (float(poly) >= 0.80):
							if (clinVarOut == 'likely_benign') or (clinVarOut == 'benign'):
								break
							else:	
								outData = Chr+'\t'+Pos+'\t'+wt+'\t'+mu+'\t'+ENST+'\t'+ENSG+'\t'+geneName+'\t'+types+'\t'+cDesc+'\t'+pDesc+'\t'+impact+'\t'+sift+'\t'+poly+'\t'+clinVarOut
								dataOut.write(outData+'\n')						
								break

	with open(outputFile) as data:
		for i, l in enumerate(data):
			pass
	print('There are {} nonsense and missense shared variants saved in {}'.format(i+1,outputFile+'.gz'))

	cmdLine = "gzip {}".format(inputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gnzip {}".format(inputFile))
		exit(1)
	
	cmdLine = "gzip {}".format(outputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(outputFile))
		exit(1)
