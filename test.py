import subprocess
import os, sys

os.system('zcat ./media/glyma.Wm82.gnm1.div.BRR2.SNPs.parents.vcf.gz | grep "#CHROM" > ./media/test.txt')