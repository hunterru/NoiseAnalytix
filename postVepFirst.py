#!/usr/local/Python/3.6.3/bin/python3

#def postVepParser(inputFile,outputFile):

import subprocess
import argparse
import os, sys, re

inputFile1 = 'OutputData/vepFolder/Tumor1_P1.vep.txt'
inputFile2 = inputFile1
outputFile = 'tempPostVepOut.tsv'

with open(inputFile1,'r') as vepVars1, open(outputFile,'w') as dataOut:
	cChr1 = '0'
	cPos1 = '0'
	cwt1 = '0'
	cmu1 = '0'
	for line1 in vepVars1:
		line1 = line1.rstrip()
		if line1[0] != "#":
			line1Sep = line1.split()
			pos1sep = line1Sep[0].split('_')
			Chr1 = pos1sep[0]
			Pos1 = pos1sep[1]
			var1 = pos1sep[2].split('/')
			wt1 = var1[0]
			mu1 = var1[1]
			if (cChr1 == Chr1) and (cPos1 == Pos1) and (cwt1 == wt1) and (cmu1 == mu1):
				continue
			else:
				with open(inputFile2,'r') as vepVars2:
					accum = []
					for line2 in vepVars2:
						line2 = line2.rstrip()
						if line2[0] != '#':
							line2Sep = line2.split()
							pos2sep = line2Sep[0].split('_')
							Chr2 = pos2sep[0]
							Pos2 = pos2sep[1]
							if (Chr1 == Chr2) and (Pos1 == Pos2):
								var2 = pos2sep[2].split('/')
								wt2 = var2[0]
								mu2 = var2[1]
								if (wt1 == wt2) and (mu1 == mu2):
									accum.append(line2Sep)
					sEnst = 0
					counter = 0
					for lines in accum:
						Enst = lines[4] 
						found1 = re.search(r"ENST0+",Enst)
						cEnst = int(Enst[len(found1.group(0)):])
						if cEnst > sEnst:
							sEnst = cEnst
							out = lines
					impactSplit = out[6].split(',')
					for s in impactSplit:
						if (s == 'stop_gained') or (s == 'missense_variant'):
							print('Keeping:',s)
							break
					cChr1 = Chr1
					cPos1 = Pos1
					cwt1 = wt1
					cmu1 = mu1

