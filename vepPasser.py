#!/usr/local/Python/3.6.3/bin/python3

def vepVars(inputFile,outputFile,vepFolder):
	import subprocess
	import argparse
	import os, sys, re
	import gzip
	
	transVepIn = vepFolder+'/tempVepIn.tsv'
	transVepOut = vepFolder+'/tempVepOut.txt'
	if inputFile.endswith('.gz') == 1:
		cmdLine = "gunzip {}".format(inputFile)
		cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cmdLine_run.returncode != 0:
			print("FAILED! Unable to gunzip {}".format(inputFile))
			exit(1)
		else:
			inputFile = inputFile[:-3]

	with open(inputFile,'r') as vepIn, open(outputFile,'w') as vepOut:
		counter = 0
		vepList = []
		for lines in vepIn:
			counter += 1
			lines = lines.rstrip()
			vepList.append(lines)
			if counter == 500:
				with open(transVepIn,'w') as tempVepIn:
					for veps in vepList:
						tempVepIn.write(veps+'\n')
				with open(transVepIn,'r') as tempVepIn:
					for i, l in enumerate(tempVepIn):
						pass
				print('In a looping block of {} input lines for VEP analysis'.format(i+1))
				vepCmd = 'vep --database -a GRCh37 --port 3337 --sift b --polyphen b --humdiv --check_existing --symbol --force -i {} -o {}'.format(transVepIn,transVepOut)
				vepCmd_run = subprocess.run(vepCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				if vepCmd_run.returncode != 0:
					print("FAILED! Unable to run VEP for variant classification.")
					exit(1)

				with open(transVepOut,'r') as vepOutput:
					for vepLines in vepOutput:
						vepLines = vepLines.rstrip()
						if vepLines[0] != '#':
							vepOut.write(vepLines+'\n')

				vepList = []
				counter = 0

		if len(vepList) > 0:			
			with open(transVepIn,'w') as tempVepIn:
				for veps in vepList:
					tempVepIn.write(veps+'\n')
			with open(transVepIn,'r') as tempVepIn:
				for i, l in enumerate(tempVepIn):
					pass
			print('Entering final block containing {} input lines for VEP analysis'.format(i+1))
			vepCmd = 'vep --database -a GRCh37 --port 3337 --sift b --polyphen b --humdiv --check_existing --symbol --force -i {} -o {}'.format(transVepIn,transVepOut)
			vepCmd_run = subprocess.run(vepCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			if vepCmd_run.returncode != 0:
				print("FAILED! Unable to run VEP for variant classification.")
				exit(1)

			with open(transVepOut,'r') as vepOutput:
				for vepLines in vepOutput:
					vepLines = vepLines.rstrip()
					vepOut.write(vepLines+'\n')

	rmCmd = 'rm {} {}'.format(transVepIn, transVepOut)
	rmCmd_run = subprocess.run(rmCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if rmCmd_run.returncode != 0:
		print("FAILED! Unable to delete temporary files during VEP passing function.")
		exit(1)

	with open(outputFile,'r') as vepFinal:
		for i, l in enumerate(vepFinal):
			pass
	print("VEP analysis complete. There are {} lines with variants for further analysis saved at {}".format(i+1, outputFile+'.gz'))

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


