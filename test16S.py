## AUTHOR: Mariana P. <m.dornelles19@gmail.com>
## USAGE: python test16S.py > 16S.fasta 
## python extract_gene_rpoB.py
## WHAT FOR? This script retrieves 16S nucleotide sequence 
## from all the .gbff files in the folder
## IMPORTANT: the script must be in the same directory of your 
## your genomes of interest
#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


files = os.listdir()
#print(files)
for f in files:
	if(f.endswith('.gbff')):
		gbank = SeqIO.parse(open(f, "r"), "genbank")
		rRNAs = []
		for genome in gbank:
			for feature in genome.features:
				if(feature.type == "rRNA"):
					product = feature.qualifiers['product'][0]
					if(product == '16S ribosomal RNA'):
						seq = feature.extract(genome.seq)
						print(">"+f.replace(".gbff",""))
						print(seq)
					break

