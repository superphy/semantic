__author__ = 'Stephen Kan'

from Bio import Entrez

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

def return_metadata_nuccore_efetch(nuccore_id):
    handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="gb", retmode="text", seq_start="1", seq_stop="1")
    return handle