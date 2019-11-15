#!/usr/local/Python/3.6.3/bin/python3

def postVepParser(inputFile,outputFile):

	import subprocess
	import argparse
	import os, sys, re

	with open(inputFile) as vepVars:
		outTemp = []
		for lines in vepVars:
			lines = lines.rstrip()
			if lines[0] != "#":
				linesSep = lines.split()
				if linesSep[6] == 'missense_variant':
					if linesSep[4] == 'ENST00000557334':
						pos0sep = linesSep[0].split('_')
						Chr = pos0sep[0]
						Pos = pos0sep[1]
						var = pos0sep[2].split('/')
						wt = var[0]
						mu = var[1]
						cDesc = 'c.'+linesSep[7]+wt+'>'+mu
						aa = linesSep[10].split('/')
						pDesc = 'p.'+aa[0]+linesSep[9]+aa[1]
						exParse = linesSep[13].split(';')
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
						out = Chr+'\t'+Pos+'\t'+wt+'\t'+mu+'\t'+geneName+'\t'+cDesc+'\t'+pDesc+'\t'+impact+'\t'+sift+'\t'+poly+'\t'+clinVarOut
						outTemp.append(out)
	outSet = set(outTemp)
	with open(outputFile,'w') as dataOut:
		dataOut.write('Chr\tPos\tWT\tMUT\tGENE\tcDesc\tpDesc\tIMPACT\tSIFT\tPoly\tClinVar\n')
		for lines in outSet:
			dataOut.write(lines+'\n')
