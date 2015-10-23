__author__ = 'Stephen Kan'

import ijson.backends.yajl2 as ijson
import sys
import traceback
import ontology_uploader
import os
import inspect

from rdflib import Graph
from eutils import return_elink_uid, return_nuccore_efetch, return_esearch_uid, only_digits
from classes import PendingGenome, generate_output
from sparql import check_genome

reload(sys)
sys.setdefaultencoding("utf-8")

genome_params = {"isolation_date":"date", "isolation_location":"location", "isolation_host":"host",
                 "isolation_source":"source"}
g = Graph()

def load_minerJSON(filename, organism):
    progress = 0
    error = 0
    dict = {}

    path = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), filename)

    with open(path, "r") as fd:
        parser = ijson.parse(fd)

        for prefix, event, value in parser:

            if ("." not in prefix) and (prefix is not "") and (event == "start_map"):
                if dict:
                    progress += 1
                    print str(progress) + ": downloading files"
                    try:
                        create_pending_genome(dict)
                    except Exception as e:
                        error = error_logging(dict, error)

                    dict.clear()

                dict.setdefault("accession", set())
                dict["name"] = prefix
                dict["accession"].add(prefix)
                dict["organism"] = organism

            if prefix.endswith(".displayname"):
                cat = prefix.split(".", 3)[1]
                dict.setdefault(cat, set())
                dict[cat].add(value)

        if dict:
            progress += 1
            print str(progress) + ": downloading files"
        try:
            create_pending_genome(dict)
        except Exception as e:
            error_logging(dict, error)


def error_logging(dict, error):
    f = open(os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), "outputs/errors.txt"), "a")
    f.write(traceback.format_exc() + "\n")
    f.write(dict["name"] + "\n" + "=======================" + "\n")
    error += 1
    print "Error %d occurred." % error
    return error


def create_pending_genome(dict):
    n = dict["name"]
    kwargs = {}

    if check_genome(n):
        print n + " already in Blazegraph."
    else:
        nuccore_id = return_esearch_uid("nuccore", n)
        (biosample, bioproject) = return_bio_ids(nuccore_id)
        kwargs.update({"bioproject":bioproject, "biosample":biosample})

        for key, value in dict.iteritems():
            if key in genome_params:
                param = genome_params.get(key)
                kwargs.update({param:value})

            elif key == "serotype":
                (Otype, Htype) = return_serotypes(value)
                kwargs.update({"Otype":Otype, "Htype":Htype})

            else:
                pass
                kwargs.update({key:value})

        PendingGenome(g, **kwargs).rdf()
        output = generate_output(g)
        ontology_uploader.upload_data(output)

def return_serotypes(serotypes):
    Otype = None
    Htype = None

    for serotype in serotypes:
        if "ONT" in serotype:
            Otype = None
        else:
            Otype = serotype.split(":")[0][1:]

        if "NM" in serotype:
            Htype = "-"
        elif "NA" in serotype:
            Htype = None
        else:
            Htype = serotype.split(":")[1][1:]

    return (Otype, Htype)


def return_bio_ids(nuccore_id):
    bioproject = set()
    biosample = set()

    try:
        bioproject = bioproject  | return_elink_uid("nuccore","bioproject",nuccore_id)
        try:
            biosample = biosample | return_elink_uid("nuccore","biosample",nuccore_id)

        except IndexError:
            for id in bioproject:
                biosample = biosample | return_elink_uid("bioproject","biosample",id)

    except IndexError:
        for record in return_nuccore_efetch(nuccore_id):
            for xref in record["GBSeq_xrefs"]:
                if xref["GBXref_dbname"] == "BioProject":
                    bioproject.add(only_digits(xref["GBXref_id"]).lstrip("0"))
                elif xref["GBXref_dbname"] == "BioSample":
                    biosample.add(only_digits(xref["GBXref_id"]).lstrip("0"))

    return (bioproject, biosample)




load_minerJSON("samples/small_pipe.json", "ecoli")
#load_minerJSON("samples/meta_pipe_result.json", "ecoli")