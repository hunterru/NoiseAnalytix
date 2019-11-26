#!/usr/local/Python/3.6.3/bin/python3

#bcftools mpileup --count-orphans --no-BAQ --regions-file ~/ErrorClass/tempBed.bed --max-depth 100000 --max-idepth 100000 --gap-frac 0.001 --per-sample-mF --ignore-RG --min-MQ 13 --fasta-ref /home/Genomes/Human/B37/human_g1k_v37_decoy_phiXAdaptr.fasta --ff UNMAP,SECONDARY,QCFAIL --annotate FORMAT/AD -Ou ~/ErrorClass/bamFolder/Tumor1_P1.bam ~/ErrorClass/bamFolder/Buffy_P1.bam | bcftools norm --fasta-ref /home/Genomes/Human/B37/human_g1k_v37_decoy_phiXAdaptr.fasta --multiallelics -any --threads 5 | grep -w -v '<\*>' | gzip > OutputData/NewOutput1.vcf.gz

#java -jar ~/BioApps/USeq_9.2.4.beta/Apps/SimpleSomaticCaller -v OutputData/NewOutput1.vcf -f 1 -t 0.00001 -x 0.6 -n 0.1 -u 1000 -a 1 -o 1000 -d 0.0001 -r 0.00001


## Example: ./aTenuWait.py myBams.txt bamFolder/ /home/Genomes/Human/B37/human_g1k_v37_decoy_phiXAdaptr.fasta ~/ErrorClass/OutputData -l 12:25398280-25398292

import subprocess
import argparse
import os, sys, re
import operator
from postSSC import postSSCparser
from fastPostParser import commonVarsNew
from vepPasser import vepVars
from postVepFirstCmd import postVepParserFirst
from postVepSecond import postVepParserSecond
from addVaf2Vep import postVepAddVaf
from varClassifier import postVepVafClass

def main():
	parser = argparse.ArgumentParser(description="""
		The 'aTenuWaitToo' program takes the general idea from aTenuWait and uses David
		Nix's Simple Somatic Caller to generate a .vcf. Freebayes seemed to work well, 
		but the time was too long for the depth being used. Also, Freebayes does not 
		account for duplicates. As with the SSC, aTenuWaitToo uses a normal reference
		(i.e., buffy coat DNA), but also uses merges potential variants from two files.
		The program the diverges and goes into VEP to assist with variant ajudication. 
		The gene, principle varient effect, and predictive effects on protein function 
		(via PolyPhen and Sift in VEP) are returned. In addition, read depth and VAF are
		also returned in a final .vcf file. Either -l or -e flag must be used (but not both!).
		""")
	parser.add_argument('bamList', help='A .txt file with a list of .bam files (and location if not in working directory). The order should be Normal, Tumor1, Tumor2')
	parser.add_argument('bamFiles', help='A directory location with the .bam files identified in bamList')
	parser.add_argument('refGenome', help='Identify reference genome (and location if not in working directory')
	parser.add_argument('bedFile', help='Target location(s) of interest formatted as a bedfile')
	parser.add_argument('output', help='Target directory for output files created by aTenuWait') 
	args = parser.parse_args()

	bamList = args.bamList
	bamFiles = args.bamFiles
	refGenome = args.refGenome
	dirOut = args.output
	bedFile = args.bedFile

	currPath = os.getcwd()
	startDir = currPath + '/'
	if os.path.exists(dirOut) == 0:
		os.mkdir(dirOut)
	mpOut = dirOut + '/mpFolder'
	if os.path.exists(mpOut) == 0:
		os.mkdir(mpOut)
	compOut = dirOut + '/compFolder'
	if os.path.exists(compOut) == 0:
		os.mkdir(compOut)
	vepOut = dirOut + '/vepFolder'
	if os.path.exists(vepOut) == 0:
		os.mkdir(vepOut)
	finalOut = dirOut + '/finalVarFolder'
	if os.path.exists(finalOut) == 0:
		os.mkdir(finalOut)

	## Taking each bam file in the user entered bamList (2 at a time) and identifying varients
	with open(bamList,'r') as bamFileList:
		bamCount = 0
		for bams in bamFileList:
	## Step 1. pass bam files to bcftools mpileup
			bamCount += 1
			print('Bam Counter at:',bamCount)			
			bams = bams.rstrip()
			if bamCount%3 == 2:
				bam1name = bams
			bamsIn = bamFiles+bams
			if bamCount%3 == 1:
				NORMAL = bamsIn
				continue
			elif bamCount%3 != 1:
				TUMOR = bamsIn
				varOut = mpOut+'/'+bams[:-4]+'.mp.vcf.gz'
				cmdLine = "bcftools mpileup --count-orphans --no-BAQ --regions-file {} --max-depth 100000 --max-idepth 100000 --gap-frac 0.001 --per-sample-mF --ignore-RG --min-MQ 13 --fasta-ref {} --ff UNMAP,SECONDARY,QCFAIL --annotate FORMAT/AD -Ou {} {} | bcftools norm --fasta-ref {} --multiallelics -any --threads 5 | grep -w -v '<\*>' | gzip > {}".format(bedFile,refGenome,TUMOR,NORMAL,refGenome,varOut)
				cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				if cmdLine_run.returncode != 0:
					print("FAILED! Unable to run bcftools mpileup for initial variant detection.")
					exit(1)	
				else:
					print("bcftools mpileup .vcf file generation complete and is located at: {}".format(varOut))
	## Step 2. pass .vcf.gz file to SimpleSomatic Caller. Example:  #java -jar ~/BioApps/USeq_9.2.4.beta/Apps/SimpleSomaticCaller -v OutputData/NewOutput1.vcf -f 1 -t 0.00001 -x 0.6 -n 0.1 -u 1000 -a 1 -o 1000 -d 0.0001 -r 5
				cmdLine = "java -jar ~/BioApps/USeq_9.2.4.beta/Apps/SimpleSomaticCaller -v {} -f 1 -t 0.00001 -x 0.6 -n 0.1 -u 1000 -a 1 -o 1000 -d 0.0001 -r 10".format(varOut)
				cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				if cmdLine_run.returncode != 0:
					print("FAILED! Unable to run SimpleSomaticCaller.")
					exit(1)
				else:
					print("SimpleSomaticCaller worked!")

				if bamCount%3 == 2:
					VAR1 = varOut[:-6]+'ssc.vcf.gz'
					print('SSC output saved at {}'.format(VAR1))
					parsVAR1 = VAR1[:-2]+'parsed.tsv'
					postSSCparser(VAR1,parsVAR1)
					CommonOut = compOut+'/'+bams[:-4]+'.shared.tsv'
					VepIn = compOut+'/'+bams[:-4]+'.vepIn.tsv'
					VepOutput = vepOut+'/'+bams[:-4]+'.vep.txt' 
					print('Parsed SSC saved at {}'.format(parsVAR1+'.gz'))
				elif bamCount%3 == 0:
					VAR2 = varOut[:-6]+'ssc.vcf.gz'
					print('SSC output saved at {}'.format(VAR2))
					parsVAR2 = VAR2[:-2]+'parsed.tsv'
					postSSCparser(VAR2,parsVAR2)
					print('Parsed SSC saved at {}'.format(parsVAR2+'.gz'))
	## Step 3. Using both parsed .vcf files, identify variants that are only present in both.
					print('Identifying common variants.')
					commonVarsNew(parsVAR1+'.gz',parsVAR2+'.gz',compOut,VepIn)

				## vepOpenCmd = 'module load vep'
				## vepOpenCmd_run = subprocess.run(vepOpenCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				## include --humdiv : HumDiv-trained model should be used for evaluating rare alleles at loci potentially involved in 
				## complex phenotypes, dense mapping of regions identified by genome-wide association studies, and analysis of natural 
				## selection from sequence data, where even mildly deleterious alleles must be treated as damaging.
				
					cmdLine = 'gunzip {}'.format(VepIn+'.gz')
					cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					if cmdLine_run.returncode != 0:
						print("FAILED! Unable to gunzip {}".format(VepIn+'.gz'))
						exit(1)
					if os.stat(VepIn).st_size > 0:
						print('Accessing VEP for potential variants. A while this may take. Patient, you must be.')
						vepVars(VepIn,VepOutput,vepOut)

						VepParsedOne = VepOutput[:-4]+'.parsed1.tsv'
						VepParsedTwo = VepOutput[:-4]+'.parsed2.tsv'
						postVepParserFirst(VepOutput+'.gz',VepParsedOne,vepOut)
						postVepParserSecond(VepParsedOne+'.gz',VepParsedTwo)

						
						adjNoise = finalOut+'/'+bam1name+'.adj_var.tsv'
						postVepAddVaf(parsVAR1+'.gz',VepParsedTwo+'.gz',adjNoise)

					elif os.stat(sharedSnvOut).st_size == 0:
						print('There are no shared variants between the two .bam files. VEP was not run on data!')
	
		sys.exit(0)

if __name__=='__main__':
	main()
