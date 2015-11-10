__author__ = 'Stephen Kan'

from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(query="samples/sampleWGS.fasta", db="../../../../Documents/BLASTDB/ValidationDB", evalue = 0.001, outfmt=5, out="results.xml")
print blastn_cline()
