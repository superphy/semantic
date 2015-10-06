__author__ = 'ubiquitin'

import json
from Bio import Entrez

f = open("parseResult.txt", "w")

def main():
    with open("samples/meta_pipe_result.json") as json_file:
        json_data = json.load(json_file)

        for accession in json_data:
            f.write("accession number: " + accession + "\n")

            Entrez.email = "stebokan@gmail.com"
            handle = Entrez.esearch(db="biosample", retmax=5, term=accession)
            record = Entrez.read(handle)
            handle.close
            f.write("Number of Biosample Files: " + record["Count"] + "\n")

            for item in record["IdList"]:
                    f.write(item + "\n")

            handle = Entrez.esearch(db="bioproject", retmax=5, term=accession)
            record = Entrez.read(handle)
            handle.close
            f.write("Number of Bioproject Files: " + record["Count"] + "\n")

            for item in record["IdList"]:
                 f.write(item + "\n")

            contents = json_data[accession]

            for item in contents:
                data = contents[item]

                for value in data:
                    try:
                        f.write(item + ": " + value["displayname"] + "\n")
                    except KeyError:

                        for some in value:
                            f.write(item + ": " + value[some]["displayname"] + "\n")
        f.write("=======================" + "\n")


main()