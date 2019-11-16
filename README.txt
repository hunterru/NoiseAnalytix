# NoiseAnalytix

usage: aTenuWait.py [-h] [-l LOCATIONS] [-e BEDFILE] bamFiles refGenome output

The 'aTenuWait' program makes an attempt to adjudicate real 
variants from noise. The program uses .bam files (aligned, 
collapsed reads) and identifies variants (SNVs, INDELs) at 
all positions via freebayes. If DNA has been sequenced in 
duplicate, then each position is compared and only variants 
present in both duplicates are preserved. Variants are then
passed to VEP for variant classification. The gene, principle 
varient effect, and predictive effects on protein function 
(via PolyPhen and Sift in VEP) are returned. In addition, read
depth and VAF are also returned in a final .vcf file. Either
-l or -e flag must be used (but not both!).

positional arguments:
  bamFiles              A .txt file with a list of .bam files (and location if
                        not in working directory
  refGenome             Identify reference genome (and location if not in
                        working directory
  output                Target directory for output files created by aTenuWait

optional arguments:
  -h, --help            show this help message and exit
  -l LOCATIONS, --locations LOCATIONS
                        Target location(s) of interest formatted as
                        chr:pos1-pos2 (ex: 12:25398280-25398292)
  -e BEDFILE, --bedfile BEDFILE
                        Target locations of interest in a bedfile

Example:
./aTenuWait.py myBams.txt /home/Genomes/Human/B37/human_g1k_v37_decoy_phiXAdaptr.fasta ~/ErrorClass/OutputData/ -e tempBed.bed


Dependents:
	postVep.py
	postVepVaf.py
	varClassifier.py[key1],key2,newPosDictSorted[key2])
