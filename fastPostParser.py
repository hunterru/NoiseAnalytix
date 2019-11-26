#!/usr/local/Python/3.6.3/bin/python3

def commonVarsNew(inputFile1,inputFile2,workDir,outputFile):

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

	tempOut1 = workDir+'/tempAwkOut1.tsv'
	cmdLine =  "awk "+"'{"+"print "+"$1,$2,$4,$5}' "+"{} > {}".format(inputFile1,tempOut1)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to pull columns from {}".format(outputFile1))
		exit(1)

	tempOut2 = workDir+'/tempAwkOut2.tsv'
	cmdLine =  "awk "+"'{"+"print "+"$1,$2,$4,$5}' "+"{} > {}".format(inputFile1,tempOut2)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to pull columns from {}".format(outputFile2))
		exit(1)

	tempOut3 = workDir+'/tempCommon.tsv'
	cmdLine = "awk "+"'NR==FNR{"+"arr[$0];next} "+"$0 in arr' {} {} > {}".format(tempOut1, tempOut2, tempOut3)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to find common variants between {} and {}".format(inputFile1, inputFile2))
		exit(1)
	else:
		count = 0
		with open(tempOut3,'r') as Data, open(outputFile,'w') as DataOut:
			for lines in Data:
				lines = lines.rstrip()
				if lines[0] != '#':
					count += 1
					lineSep = lines.split()
					vepIn = lineSep[0]+'\t'+lineSep[1]+'\t.\t'+lineSep[2]+'\t'+lineSep[3]+'\t.\tPASS\t.\n'
					DataOut.write(vepIn)
		print('There were {} common variants identified and saved in {}'.format(count,outputFile))

	cmdLine = "rm {} {} {}".format(tempOut1, tempOut2, tempOut3)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to remove temporary files {} and {}".format(tempOut1,tempOut2))
		exit(1)

	cmdLine = "gzip {}".format(outputFile)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(outputFile))
		exit(1)

	cmdLine = "gzip {}".format(inputFile1)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(inputFile1))
		exit(1)
	cmdLine = "gzip {}".format(inputFile2)
	cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if cmdLine_run.returncode != 0:
		print("FAILED! Unable to gzip {}".format(inputFile2))
		exit(1)
