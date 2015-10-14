__author__ = 'Stephen Kan'

import json
import string
import sys
import traceback

from Bio import Entrez
from superphy_classes import PendingGenome, generate_output

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
Otype = None
Htype = None

metadata = [nuccore, bioproject, biosample, syndrome, strain, serotype, isolation_location, isolation_date, isolation_host, isolation_source]


def main():
    with open("samples/small_pipe.json") as json_file:
        json_data = json.load(json_file)
        i = 0
        j = 0

        for accession in json_data:
            for metadatum in metadata:
                metadatum.clear()
            global Otype, Htype
            Otype = None
            Htype = None

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

            if serotype is not None:
                for sero in serotype:
                    if "ONT" in sero:
                        global Otype
                        Otype = None
                    else:
                        global Otype
                        Otype = sero.split(":")[0][1:]

                    if "NM" in sero:
                        global Htype
                        Htype = "H-"
                    elif "NA" in sero:
                        global Htype
                        Htype = "H-"
                    else:
                        global Htype
                        Htype = sero.split(":")[1][1:]

            try:
                kwargs = {"date": isolation_date, "location": isolation_location, "accession": nuccore,
                          "bioproject": bioproject, "biosample": biosample, "strain": strain, "organism": "ecoli",
                          "from_host": isolation_host, "from_source": isolation_source, "syndrome": syndrome,
                          "Otype": Otype, "Htype": Htype}
                PendingGenome(name, **kwargs).rdf()
            except Exception as e:
                f = open("outputs/errors.txt", "a")
                f.write(traceback.format_exc() + "\n")
                f.write(accession + "\n" + "=======================" + "\n")
                j+=1
                print "Error %d occurred." %j

    generate_output("outputs/results.ttl")
    print "%d errors occurred" %j




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