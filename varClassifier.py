#!/usr/local/Python/3.6.3/bin/python3

def postVepAddVaf(countFile,vepFile,outputFile):

	import subprocess
	import argparse
	import os, sys, re

	#countFile = 'OutputData/Tumor1.counts.tsv'
	#vepFile = 'OutputData/Tumor1.sharedSnvVars.vepOutput.parsed.txt'

	with open(countFile,'r') as COUNTS, open(outputFile,'w') as DATAOUT:
		DATAOUT.write('Chr\tPos\tWT\tMUT\tGENE\tcDesc\tpDesc\tDP\tvarCount\tpVAF\tIMPACT\tSIFT\tPoly\tClinVar\n')
		for countlines in COUNTS:
			countlines = countlines.rstrip()
			countSplit = countlines.split('\t')
			if countSplit[0] != 'Chr':
				with open(vepFile,'r') as VART:
					for varlines in VART:
						varlines = varlines.rstrip()
						varSplit = varlines.split('\t')
						if varSplit[0] != 'Chr':
							if countSplit[0] == varSplit[0]:
								if countSplit[1] == varSplit[1]:
									dp = countSplit[2]
									if varSplit[3] == 'C':
										varDp = countSplit[3]
									elif varSplit[3] == 'A':
										varDp = countSplit[4]
									elif varSplit[3] == 'T':
										varDp = countSplit[5]
									elif varSplit[3] == 'G':
										varDp = countSplit[6] 
									out = varSplit[0:7]
									out.append(dp)
									out.append(varDp)
									vaf = float(varDp) / float(dp) * 100
									out.append('{0:.{1}f}'.format(vaf,3))
									out.append(varSplit[7])
									out.append(varSplit[8])
									out.append(varSplit[9])
									if len(varSplit) == 11:
										out.append(varSplit[10])
										DATAOUT.write(out[0]+'\t'+out[1]+'\t'+out[2]+'\t'+out[3]+'\t'+out[4]+'\t'+out[5]+'\t'+out[6]+'\t'+out[7]+'\t'+out[8]+'\t'+out[9]+'\t'+out[10]+'\t'+out[11]+'\t'+out[12]+'\t'+out[13]+'\n')
									else:
										DATAOUT.write(out[0]+'\t'+out[1]+'\t'+out[2]+'\t'+out[3]+'\t'+out[4]+'\t'+out[5]+'\t'+out[6]+'\t'+out[7]+'\t'+out[8]+'\t'+out[9]+'\t'+out[10]+'\t'+out[11]+'\t'+out[12]+'\n')
