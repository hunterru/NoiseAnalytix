# NoiseAnalytix

The 'aTenuWait' program makes an attempt to adjudicate real variants from
noise. The program uses .bam files (aligned, collapsed reads) and identifies
variants (SNVs, INDELs) at all positions via freebayes. If DNA has been
sequenced in duplicate, then each position is compared and only variants
present in both duplicates are preserved. Variants are then passed to VEP for
variant classification. The gene, principle varient effect, and predictive
effects on protein function (via PolyPhen and Sift in VEP) are returned. In
addition, read depth and VAF are also returned in a final .vcf file.

optional arguments:
  -h, --help            show this help message and exit
  -b BAMFILES, --bamFiles BAMFILES
                        A .txt file with a list of .bam files (and location if
                        not in working directory
  -d {True,False}, --duplicates {True,False}
                        Was your sequencing data done in duplicate?
  -f REFGENOME, --refGenome REFGENOME
                        Identify reference genome (and location if not in
                        working directory
  -l LOCATIONS, --locations LOCATIONS
                        Target location(s) of interest or provide a bedfile
  -o OUTPUT, --output OUTPUT
                        Target directory for output files created by aTenuWait

Example of use:
./aTenuWait.py -b myBams.txt -d True -f /home/Genomes/Human/B37/human_g1k_v37_decoy_phiXAdaptr.fasta -l 12:25398280-25398292 -o ~/ErrorClass/OutputData/

Dependent files:
postVep.py
postVepVaf.py
varClassifier.py
