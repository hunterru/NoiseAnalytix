#!/usr/local/Python/3.6.3/bin/python3

## Example: ./aTenuWait.py -b myBams.txt -d True -f /home/Genomes/Human/B37/human_g1k_v37_decoy_phiXAdaptr.fasta -l 12:25398280-25398292 -o ~/ErrorClass/OutputData/

import subprocess
import argparse
import os, sys, re
import operator
from postVep import postVepParser
from postVepVaf import postVepAddVaf

def main():
	parser = argparse.ArgumentParser(description="""
		The 'aTenuWait' program makes an attempt to adjudicate real variants from noise. 
		The program uses .bam files (aligned, collapsed reads) and identifies variants 
		(SNVs, INDELs) at all positions via freebayes. If DNA has been sequenced in 
		duplicate, then each position is compared and only variants present in both
		duplicates are preserved. Variants are then passed to VEP for variant classification.  
		The gene, principle varient effect, and predictive effects on protein function 
		(via PolyPhen and Sift in VEP) are returned. In addition, read depth and VAF are
		also returned in a final .vcf file.
		""")
	parser.add_argument('-b', '--bamFiles', required=True, help='A .txt file with a list of .bam files (and location if not in working directory')
	parser.add_argument('-d', '--duplicates', choices=['True','False'], required=True, help='Was your sequencing data done in duplicate?')
	parser.add_argument('-f', '--refGenome', required=True, help='Identify reference genome (and location if not in working directory') 
	parser.add_argument('-l', '--locations', required=True, help='Target location(s) of interest or provide a bedfile')
	parser.add_argument('-o', '--output', required=True, help='Target directory for output files created by aTenuWait')
	args = parser.parse_args()

	bamFiles = args.bamFiles
	refGenome = args.refGenome
	locations = args.locations
	dirOut = args.output

	with open(bamFiles,'r') as bamList:
		bamCount = 0
		for bams in bamList:
			bams = bams.rstrip()
			bamCount += 1
			varOut = bams[:-4]+'.fb.var.vcf'
			cmdLine = "freebayes -f {} -r {} -C 1 -F 0 -K -X -u --report-monomorphic {} > {}".format(refGenome, locations, bams, dirOut+varOut)
			cmdLine_run = subprocess.run(cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			if cmdLine_run.returncode != 0:
				print("FAILED! Unable to run Freebayes for variant detection.")
				exit(1)	
			else:
				print("freebayes generation of {} complete!\n".format(varOut))

			dataIn = dirOut+varOut
			snvVarOut = dirOut+varOut[:-4]+'.pos.snv.vcf'				#fileName of SNVs to be checked against the other duplicate
			indelVarOut = dirOut+varOut[:-4]+'.pos.indel.vcf'		#fileName of INDELs to be checked against the other duplicate
			totalCounts = dirOut+bams[:-4]+'.counts.tsv'		#fileName for data containing read depth information at each position for each base
			with open(dataIn,'r') as VARS, open(snvVarOut,'w') as snvPos, open(indelVarOut,'w') as indelPos, open(totalCounts,'w') as countsPos:
				countsPos.write('Chr\tPos\tTotal\tC\tA\tT\tG\n')
				for lines in VARS:
					lines = lines.rstrip()
					if lines[0] != "#":
						info = lines.split()
						Chr = info[0]
						Pos = info[1]
						wt = info[3]
						alts = info[4].split(',')
						if alts[0] == '.':												#the case where an 'alt' is absent (only WT)
							varSearch = info[7].split(';')
							for data in varSearch:
								if data.startswith('DP='):
									dp = data[3:]
									if wt == 'C':
										countsPos.write('{}\t{}\t{}\t{}\t0\t0\t0\n'.format(Chr,Pos,dp,dp))
									elif wt == 'A':
										countsPos.write('{}\t{}\t{}\t0\t{}\t0\t0\n'.format(Chr,Pos,dp,dp))
									elif wt == 'T':
										countsPos.write('{}\t{}\t{}\t0\t0\t{}\t0\n'.format(Chr,Pos,dp,dp))
									elif wt == 'G':
										countsPos.write('{}\t{}\t{}\t0\t0\t0\t{}\n'.format(Chr,Pos,dp,dp))
						else:																			#an 'alt' is present
							varSearch = info[7].split(';')
							for data in varSearch:
								if data.startswith('DP='):
									dp = data[3:]
								if data.startswith('RO='):
									dpWT = data[3:]
								if data.startswith('AO='):
									temp = data[3:]
									dpAlt = temp.split(',')	
								if data.startswith('TYPE='):
									temp = data[5:]
									typeAlt = temp.split(',')

							ntCounts = 0
							baseDict = dict()
							for nt in wt:
								newPos = str(int(Pos) + ntCounts)
								if nt == 'C':
									baseDict[newPos] = [str(dp), str(dpWT), str(0), str(0), str(0)]
								elif nt == 'A':
									baseDict[newPos] = [str(dp), str(0), str(dpWT), str(0), str(0)]
								elif nt == 'T':
									baseDict[newPos] = [str(dp), str(0), str(0), str(dpWT), str(0)]
								elif nt == 'G':
									baseDict[newPos] = [str(dp), str(0), str(0), str(0), str(dpWT)]
								ntCounts += 1

							varCounts = 0
							for var in alts:
								if typeAlt[varCounts] == 'snp':
									ntCounts = 0
									for nt in var:
										newPos = str(int(Pos) + ntCounts)
										if nt != wt[ntCounts]:
											if nt == 'C':
												baseDict[newPos][1] = dpAlt[varCounts]
											elif nt == 'A':
												baseDict[newPos][2] = dpAlt[varCounts]
											elif nt == 'T':
												baseDict[newPos][3] = dpAlt[varCounts]
											elif nt == 'G':
												baseDict[newPos][4] = dpAlt[varCounts] 
											snvPos.write(Chr+'\t'+newPos+'\t'+wt[ntCounts]+'\t'+nt+'\n')
											varCounts += 1
											break
										else:
											ntCounts += 1
								elif typeAlt[varCounts] == 'del':
									indelPos.write(Chr+'\t'+newPos+'\t'+wt+'\t'+var+'\n')
									varCounts += 1
								elif typeAlt[varCounts] == 'ins':
									indelPos.write(Chr+'\t'+newPos+'\t'+wt+'\t'+var+'\n')
									varCounts += 1

							baseDictSorted = sorted(baseDict.items(), key=operator.itemgetter(0))
							for values in baseDictSorted:
								temp = values[1]
								countsPos.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(Chr,values[0],temp[0],temp[1],temp[2],temp[3],temp[4]))	
			print('Count files and variant files written!')

			if bamCount == 1:
				compSnvFile1 = snvVarOut
				compIndelFile1 = indelVarOut
				sharedSnvOut = dirOut+bams[:-4]+'.sharedSnvVars.vcf'
				sharedIndelOut = dirOut+bams[:-4]+'.sharedIndelVars.vcf'
				bam1counts = totalCounts

			if bamCount == 2:
				print("Comparing variants from .bam files.")
				with open(compSnvFile1,'r') as snvMod1, open(sharedSnvOut,'w') as shareVar:
					for vars1 in snvMod1:
						vars1 = vars1.rstrip()
						with open(snvVarOut,'r') as snvMod2:
							for vars2 in snvMod2:
								vars2 = vars2.rstrip()
								if vars1 == vars2:
									temp = vars1.split()
									shareVar.write(temp[0]+'\t'+temp[1]+'\t'+'.'+'\t'+temp[2]+'\t'+temp[3]+'\t'+'.'+'\t'+'PASS'+'\t'+'.\n')
									#print(temp[0]+'\t'+temp[1]+'\t'+'.'+'\t'+temp[2]+'\t'+temp[3]+'\t'+'.'+'\t'+'PASS'+'\t'+'.\n')
				print('The shared variant file has been written as {}.'.format(sharedSnvOut))
	
	#vepOpenCmd = 'module load vep'
	#vepOpenCmd_run = subprocess.run(vepOpenCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# include --humdiv : HumDiv-trained model should be used for evaluating rare alleles at loci potentially involved in complex phenotypes, dense mapping of regions identified by genome-wide association studies, and analysis of natural selection from sequence data, where even mildly deleterious alleles must be treated as damaging.
	vepOutput = sharedSnvOut[:-4]+'.vepOutput.txt'
	vepCmd = 'vep --database -a GRCh37 --port 3337 --sift b --polyphen b --humdiv --check_existing --symbol --force -i {} -o {}'.format(sharedSnvOut,vepOutput)
	vepCmd_run = subprocess.run(vepCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
	if vepCmd_run.returncode != 0:
		print("FAILED! Unable to run VEP for variant classification.")
		exit(1)
	else:
		print("vepOutput.txt file generation complete!\n")

	vepParsed = vepOutput[:-4]+'.parsed.txt'
	postVepParser(vepOutput,vepParsed)
	print('The vepOutput.txt file has been written and saved as {}'.format(vepParsed))

	vepVaf = vepParsed[:-4]+'.VAF.txt'
	postVepAddVaf(bam1counts,vepParsed,vepVaf)
	print('The variant allele frequency information has been added and saved as {}'.format(vepVaf))


	sys.exit(0)


#vep --database -a GRCh37 --port 3337 --sift b --polyphen b --humdiv --check_existing --symbol --force -i sharedVars.vcf -o vepOutput.txt

if __name__=='__main__':
	main()
