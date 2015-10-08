__author__ = 'ubiquitin'

import json
import string
import sys
import object_to_rdf_converter
from Bio import Entrez

reload(sys)
sys.setdefaultencoding("utf-8")

nuccore = set()
bioproject = set()
biosample = set()
syndrome = set()
strain = set()
serotype = set()
isolation_location = set()
isolation_host = set()
isolation_date = set()
isolation_source = set()

metadata = [nuccore, bioproject, biosample, syndrome, strain, serotype, isolation_location, isolation_date, isolation_host, isolation_source]


def main():
    with open("samples/small_pipe.json") as json_file:
        json_data = json.load(json_file)
        i = 0

        for accession in json_data:
            for metadatum in metadata:
                metadatum.clear()

            i+=1
            print str(i) + ": downloading files"
            nuccore.add(accession)
            name = accession
            get_DBlink(get_nuccore_id(accession))

            contents = json_data[accession]

            for item in contents:
                data = contents[item]

                for value in data:
                    if value is not None:
                        try:
                           eval(item).add(value["displayname"])
                        except KeyError:
                            for some in value:
                                eval(item).add(value[some]["displayname"])

            object_to_rdf_converter.create_pending_genome(name=name, date=isolation_date, location=isolation_location, accession=nuccore,
                                                          bioproject=bioproject, biosample=biosample, strain=strain, organism="ecoli")

    object_to_rdf_converter.generate_output()




def get_nuccore_id(accession):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.esearch(db="nuccore", retmax=5, term=accession)
    record = Entrez.read(handle)

    for id in record["IdList"]:
        return id

def get_DBlink(nuccore_id):
    BPid = None
    Entrez.email = "stebokan@gmail.com"
    try:
        handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=nuccore_id)
        records = Entrez.parse(handle)

        for record in records:
            bioproject.add(record["LinkSetDb"][0]["Link"][0]["Id"])
            BPid = record["LinkSetDb"][0]["Link"][0]["Id"]

        try:
            handle = Entrez.elink(dbfrom="nuccore", db="biosample", id=nuccore_id)
            records = Entrez.parse(handle)

            for record in records:
                biosample.add(record["LinkSetDb"][0]["Link"][0]["Id"])

        except IndexError:
            handle = Entrez.elink(dbfrom="bioproject", db="biosample", id=BPid)
            records = Entrez.parse(handle)

            for record in records:
                biosample.add(record["LinkSetDb"][0]["Link"][0]["Id"])

    except IndexError:
        handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="gb", retmode="xml")
        records = Entrez.parse(handle)

        for record in records:
            bioproject.add(only_digits(record["GBSeq_xrefs"][0]["GBXref_id"]))
            biosample.add(only_digits(record["GBSeq_xrefs"][1]["GBXref_id"]))


def only_digits(str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.digits)
    return str.translate(all, nodigs)

main()