__author__ = 'ubiquitin'

import json
import string
from Bio import Entrez
import sys

reload(sys)
sys.setdefaultencoding("utf-8")
f = open("outputs/parseResult.txt", "w")

"""
syndrome
strain
serotype
isolation_location
isolation_host
isolation_date
isolation_source
"""


def main():
    with open("samples/meta_pipe_result.json") as json_file:
        json_data = json.load(json_file)
        i = 0
        duplication = 0

        for accession in json_data:
            f.write("accession number: " + accession + "\n")

            i+=1
            print str(i) + ": downloading files"

            get_DBlink(get_nuccore_id(accession))

            contents = json_data[accession]

            for item in contents:
                data = contents[item]

                for value in data:
                    if value is not None:
                        try:
                            f.write(item + ": " + value["displayname"] + "\n")
                        except KeyError:
                            duplication+=1
                            for some in value:
                                f.write(item + ": " + value[some]["displayname"] + "\n")


            f.write("=======================" + "\n")
    print duplication

def get_nuccore_id(accession):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.esearch(db="nuccore", retmax=5, term=accession)
    record = Entrez.read(handle)

    for id in record["IdList"]:
        return id

def get_DBlink(nuccore_id):
    BPid = None
    try:
        Entrez.email = "stebokan@gmail.com"
        handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=nuccore_id)
        records = Entrez.parse(handle)

        for record in records:
            f.write("BioProject Id: " + record["LinkSetDb"][0]["Link"][0]["Id"] + "\n")
            BPid = record["LinkSetDb"][0]["Link"][0]["Id"]

        try:
            handle = Entrez.elink(dbfrom="nuccore", db="biosample", id=nuccore_id)
            records = Entrez.parse(handle)

            for record in records:
                f.write("BioSample Id: " + record["LinkSetDb"][0]["Link"][0]["Id"] + "\n")
        except IndexError:
            handle = Entrez.elink(dbfrom="bioproject", db="biosample", id=BPid)
            records = Entrez.parse(handle)

            for record in records:
                f.write( "BioSample Id: " + record["LinkSetDb"][0]["Link"][0]["Id"] + "\n")
    except IndexError:
        Entrez.email = "stebokan@gmail.com"
        handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="gb", retmode="xml")
        records = Entrez.parse(handle)

        for record in records:
            f.write("BioProject Id: " + only_digits(record["GBSeq_xrefs"][0]["GBXref_id"]) + "\n")
            f.write("BioSample Id: " + only_digits(record["GBSeq_xrefs"][1]["GBXref_id"]) + "\n")


def only_digits(str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.digits)
    return str.translate(all, nodigs)

main()