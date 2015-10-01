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
#print(first_record.annotations['references'])

#for reference in first_record.annotations['references']:
#    print reference