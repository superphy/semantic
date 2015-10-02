__author__ = 'ubiquitin'

from Bio import SeqIO

"""
NOTE: WITH GENBANK FILES, you must select "show sequence" before
downloading to get the actual sequence for SeqRecord in Biopython
"""


record_iterator = SeqIO.parse("Sakai.gb", "genbank")
first_record = next(record_iterator)
#print(first_record)
print(first_record.annotations)
print(first_record.annotations.keys())
print(first_record.annotations['organism'])
print(first_record.annotations['source'])
#print(first_record.annotations['references'])

organism = (first_record.annotations['organism'])

#for reference in first_record.annotations['references']:
#    print reference

"""
namespace = "https://www.github.com/superphy/"
name = "O157:H7_Str._Sakai"
OType = "O157"
HType = "H7"
strain = "Sakai"
organism = "Ecoli"
isolation_date = "1996"
isolation_from_host = "from_Human"
isolation_from_source = "feces"
isolation_syndrome = "hemorrhagic_colitis"
accession = "NC_002695"
"""