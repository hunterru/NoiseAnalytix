#!/usr/local/Python/3.6.3/bin/python3

def postVepAddVaf(countFile,vepFile,outputFile):

	import subprocess
	import argparse
	import os, sys, re

	if countFile.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(countFile)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(countFile))
			exit(1)
		else:
			countFile = countFile[:-3]

	if vepFile.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(vepFile)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(vepFile))
			exit(1)
		else:
			vepFile = vepFile[:-3]

	with open(vepFile,'r') as veps, open(outputFile,'w') as dataOut:
		dataOut.write('Chr\tPos\tType\tWT\tMUT\tVAR_DP\tVAR_VAF\tWT_DP\tWT_VAF\tENST\tENSG\tGENE\tcDesc\tpDesc\tIMPACT\tSIFT\tPoly\tClinVar\n')
		for veplines in veps:
			veplines = veplines.rstrip()
			vepSplit = veplines.split()
			if vepSplit[0] != 'Chr':
				Chr = vepSplit[0]
				Pos = vepSplit[1]
				wt = vepSplit[2]
				mu = vepSplit[3]
				with open(countFile,'r') as counters:
					for countlines in counters:
						countlines = countlines.rstrip()
						countSplit = countlines.split()
						if countSplit[0] != '#Chr':
							if (countSplit[0] == Chr) and (countSplit[1] == Pos) and (countSplit[3] == wt) and (countSplit[4] == mu):
								counter = 0
								for bits in vepSplit:
									if counter == 0:
										out = bits
										counter += 1
									elif counter == 4:
										out = out+'\t'+countSplit[5]+'\t'+countSplit[6]+'\t'+countSplit[7]+'\t'+countSplit[8]
										counter += 1
									else:
										out = out+'\t'+bits
										counter += 1
								dataOut.write(out+'\n')
								break

	cmdLine = "gzip {}".format(countFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(countFile))
		exit(1)

	cmdLine = "gzip {}".format(vepFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(vepFile))
		exit(1)

	with open(outputFile) as dataOut:
		for i,l in enumerate(dataOut):
			pass

	cmdLine = "gzip {}".format(outputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(outputFile))
		exit(1)

	print('DONE! The final {} adjudicated variants with VAF has been saved at {}'.format(i+1,outputFile+'.gz'))
