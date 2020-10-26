## AUTHOR: Mariana P. <m.dornelles19@gmail.com>
## USAGE: python extract_gene_rpoB.py > rboB.fasta
## WHAT FOR? This script retrieves the nucleotide sequence 
## of the gene of interest in all .gbff files in the folder.
## HOW TO CHANGE THE GENE? Go to line 37 and replace the current
## product for your product of interest. Then, go to line 38 and 
## replace the gene name for your gene name of interest. 
## IMPORTANT: the script must be in the same directory of your 
## your genomes of interest
from Bio import SeqIO
import sys
import os

def get_gene_feature_with_qualifier_value(seq_record, name, value):
    """Function to look for CDS feature by annotation value in sequence record.
    
    e.g. You can use this for finding features by locus tag, gene ID, or protein ID.
    """
    # Loop over the features
    for feature in genome_record.features:
        if feature.type == "gene" and value in feature.qualifiers.get(name, []):
            return feature
    # Could not find it
    return None

def get_cds_feature_with_qualifier_value(seq_record, name, value):
    """Function to look for CDS feature by annotation value in sequence record.
    
    e.g. You can use this for finding features by locus tag, gene ID, or protein ID.
    """
    # Loop over the features
    for feature in genome_record.features:
        if feature.type == "CDS" and value in feature.qualifiers.get(name, []):
            return feature
    # Could not find it

    return None


##SET HERE THE PRODUCT AND GENE NAME 
##OF YOU GENE OF INTEREST
product="DNA-directed RNA polymerase subunit beta"
gene="rpoB"


files = os.listdir()
#print(files)
for f in files:
	if(f.endswith('.gbff')):
		genome=f
		genome_record = SeqIO.read(genome, "genbank")
		cds_feature = get_cds_feature_with_qualifier_value(genome_record, "product", product)

#		print("CDS")
		if cds_feature is not None:
			taxonomy = genome_record.annotations['taxonomy']
			gene_sequence = cds_feature.extract(genome_record.seq)
			protein_sequence = gene_sequence.translate(table=11, cds=True)
#			print(cds_feature)
			print(">"+genome.replace(".gbff",""))
			print(gene_sequence)
			continue
		else:
			print(">"+genome.replace(".gbff",""))
			print("does not have this cds feature")
	
		gene_feature = get_gene_feature_with_qualifier_value(genome_record, "gene", gene)
#		print("GENE")
		if gene_feature is not None: 
#			print(gene_feature)
			taxonomy = genome_record.annotations['taxonomy']
			gene_sequence = gene_feature.extract(genome_record.seq)
			protein_sequence = gene_sequence.translate(table=11, cds=True)
#			print(gene_feature)
			print(">"+genome.replace(".gbff",""))
			print(gene_sequence)
			continue
		else:
			print(">"+genome.replace(".gbff",""))
			print("does not contain this gene feature")

