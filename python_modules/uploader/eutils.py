__author__ = 'Stephen Kan'

from Bio import Entrez
import string

Entrez.email = "stebokan@gmail.com"


def return_esearch_uid(db, term):
    handle = Entrez.esearch(db=db, retmax=5, term=term)
    record = Entrez.read(handle)
    return str(record["IdList"][0])

def return_elink_uid(dbfrom, db, id):
    handle = Entrez.elink(dbfrom=dbfrom, db=db, id=id)
    records = Entrez.parse(handle)
    ret = set()

    for record in records:
        ret.add(record["LinkSetDb"][0]["Link"][0]["Id"])

    return ret

def return_nuccore_efetch(nuccore_id):
    handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="gb", retmode="xml")
    return Entrez.parse(handle)


def only_digits(str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.digits)
    return str.translate(all, nodigs)


"""
print return_esearch_uid("nuccore", "ANVW00000000")
print return_elink_uid("nuccore", "bioproject", "431333453")
print return_elink_uid("nuccore", "biosample", return_esearch_uid("nuccore", "ANVW00000000"))

for record in return_nuccore_efetch("431333453"):
    print only_digits(record["GBSeq_xrefs"][0]["GBXref_id"]).lstrip("0")
    print only_digits(record["GBSeq_xrefs"][1]["GBXref_id"]).lstrip("0")
"""