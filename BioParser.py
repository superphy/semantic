__author__ = 'ubiquitin'

import json
from Bio import Entrez
import sys

reload(sys)
sys.setdefaultencoding("utf-8")
f = open("parseResult.txt", "w")

def main():
    with open("samples/meta_pipe_result.json") as json_file:
        json_data = json.load(json_file)

        for accession in json_data:
            f.write("accession number: " + accession + "\n")
            print "hello world"

            get_DBlink(get_nuccore_id(accession))

            contents = json_data[accession]

            for item in contents:
                data = contents[item]

                for value in data:
                    if value is not None:
                        try:
                            f.write(item + ": " + value["displayname"] + "\n")
                        except KeyError:
                            for some in value:
                                f.write(item + ": " + value[some]["displayname"] + "\n")


            f.write("=======================" + "\n")


def get_nuccore_id(accession):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.esearch(db="nuccore", retmax=5, term=accession)
    record = Entrez.read(handle)

    for id in record["IdList"]:
        return id

def get_DBlink(nuccore_id):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=nuccore_id)
    records = Entrez.parse(handle)

    for record in records:
        f.write("BioProject Id: " + record["LinkSetDb"][0]["Link"][0]["Id"] + "\n")

    handle = Entrez.elink(dbfrom="nuccore", db="biosample", id=nuccore_id)
    records = Entrez.parse(handle)

    for record in records:
        f.write("BioSample Id: " + record["LinkSetDb"][0]["Link"][0]["Id"] + "\n")

main()