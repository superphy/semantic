__author__ = 'ubiquitin'

import json
from Bio import Entrez

def count_entries():

    i = 0

    with open("samples/meta_pipe_result.json") as json_file:
        json_data = json.load(json_file)

        for accession in json_data:
            i += 1

    print i

def get_nuccore_id(accession):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.esearch(db="nuccore", retmax=5, term=accession)
    record = Entrez.read(handle)

    print ("Number of Nuccore Files: " + record["Count"])

    for id in record["IdList"]:
        print (id)

def get_DBfile(nuccore_ids):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.efetch(db="nuccore", id=nuccore_ids, rettype="gb", retmode="xml")
    records = Entrez.parse(handle)

    for record in records:
        print record["GBSeq_xrefs"][0]["GBXref_id"]
        print record["GBSeq_xrefs"][1]["GBXref_id"]

def get_DBlink(nuccore_ids):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=nuccore_ids)
    records = Entrez.parse(handle)

    for record in records:
        print "BioProject Id: " + record["LinkSetDb"][0]["Link"][0]["Id"]

    handle = Entrez.elink(dbfrom="nuccore", db="biosample", id=nuccore_ids)
    records = Entrez.parse(handle)

    for record in records:
        print "BioSample Id: " + record["LinkSetDb"][0]["Link"][0]["Id"]

get_nuccore_id("CAFL00000000")
get_nuccore_id("JHIX00000000")
get_nuccore_id("CP001855")

get_DBfile("CAFL00000000")
get_DBlink("312944605")