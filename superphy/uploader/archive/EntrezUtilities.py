__author__ = 'Stephen Kan'


"""
It is probably best to do start truncating Biosample and Bioproject files so that they do not contain alphabetical characters
"""

from Bio import Entrez

terms = {"H112180283", "UCD_JA17", "2-210-07_S1_C2", "UMEA 3304-1", "E851/71"}

for term in terms:
    print "Searching for " + term

    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.esearch(db="biosample", retmax=5, term=term)
    record = Entrez.read(handle)
    handle.close
    print "Number of Biosample Files: " + record["Count"]

    for item in record["IdList"]:
        print item

    handle = Entrez.esearch(db="bioproject", retmax=5, term=term)
    record = Entrez.read(handle)
    handle.close
    print "Number of Bioproject Files: " + record["Count"]

    for item in record["IdList"]:
        print item

    handle = Entrez.esearch(db="nucleotide", retmax=5, term=term)
    record = Entrez.read(handle)
    handle.close
    print  "Number of Genbank Files: " + record["Count"]

    for item in record["IdList"]:
        print item
