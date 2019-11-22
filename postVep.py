#!/usr/local/Python/3.6.3/bin/python3

def postVepParser(inputFile,outputFile):

	import subprocess
	import argparse
	import os, sys, re

	with open(inputFile,'r') as vepVars1, open(outputFile,'w') as dataOut:
		dataOut.write('Chr\tPos\tWT\tMUT\tGENE\tcDesc\tpDesc\tIMPACT\tSIFT\tPoly\tClinVar\n')
		Chr1dup = '0'
		Pos1dup = '0'
		mu1dup = '0'
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
				Enst1 = line1Sep[4]
				found1 = re.search(r"ENST0+",Enst1)
				cENST1 = int(Enst1[len(found1.group(0)):])
				if (Chr1 == Chr1dup) and (Pos1 == Pos1dup) and (mu1 == mu1dup):
					continue
				else:
					with open(inputFile,'r') as vepVars2:
						for line2 in vepVars2:
							line2 = line2.rstrip()
							if line2[0] != "#":
								line2Sep = line2.split()
								pos2sep = line2Sep[0].split('_')
								Chr2 = pos2sep[0]
								if Chr1 == Chr2:
									Pos2 = pos2sep[1]
									if Pos1 == Pos2:
										var2 = pos2sep[2].split('/')
										wt2 = var2[0]
										mu2 = var2[1]
										if mu1 == mu2:
											Enst2 = line2Sep[4]
											found2 = re.search(r"ENST0+",Enst2)
											cENST2 = int(Enst2[len(found2.group(0)):])
											if cENST1 >= cENST2:
												linesSep = line1Sep
											elif cENST2 > cENST1:
												linesSep = line2Sep
					if (linesSep[6] == 'missense_variant') or (linesSep[6] == 'stop_gained'):
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
						if linesSep[6] == 'missense_variant':	
							outData = Chr+'\t'+Pos+'\t'+wt+'\t'+mu+'\t'+geneName+'\t'+cDesc+'\t'+pDesc+'\t'+impact+'\t'+sift+'\t'+poly+'\t'+clinVarOut
						elif linesSep[6] == 'stop_gained':
							outData = Chr+'\t'+Pos+'\t'+wt+'\t'+mu+'\t'+geneName+'\t'+cDesc+'\t'+pDesc+'\t'+impact
						dataOut.write(outData+'\n')
						Chr1dup = Chr1
						Pos1dup = Pos1
						mu1dup = mu1 
