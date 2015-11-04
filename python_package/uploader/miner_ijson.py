__author__ = 'Stephen Kan'

from ijson.backends import YAJLImportError
try:
    import ijson.backends.yajl2 as ijson
except YAJLImportError:
    import ijson.backends.yajl as ijson

import sys
import traceback
import os
import inspect

from rdflib import Graph
from eutils import return_elink_uid, return_esearch_uid
from classes import PendingGenome, generate_output
from sparql import check_NamedIndividual
from ontology_uploader import upload_data

reload(sys)
sys.setdefaultencoding("utf-8")


class MinerDataUploader(object):
    genome_params = {"isolation_date":"date", "isolation_location":"location", "isolation_host":"host",
                     "isolation_source":"source"}
    currdir = os.path.dirname(inspect.getfile(inspect.currentframe()))

    def __init__(self, filename, organism):
        self.progress = 0
        self.error = 0
        self.filename = filename
        self.organism = organism
        self.dict = {}


    def upload(self):
        self.load_JSON()
        self.iterate()
        print "%d genomes parsed, %d errors occurred." %(self.progress, self.error)


    def load_JSON(self):
        path = os.path.join(self.currdir, self.filename)
        fd = open(path, "r")
        self.data = ijson.parse(fd)


    def iterate(self):
        for prefix, event, value in self.data:
            if ("." not in prefix) and (prefix is not "") and (event == "start_map"):
                self.add_to_graph()
                self.start_new_genome(prefix)

            if prefix.endswith(".displayname"):
                self.add_genome_parameter(prefix, value)

        self.add_to_graph()


    def add_to_graph(self):
        if self.dict:
            try:
                self.progress += 1
                print str(self.progress) + ": downloading files"
                self.create_pending_genome()
            except Exception as e:
                self.error_logging()

            self.dict.clear()


    def start_new_genome(self, prefix):
        self.dict.setdefault("accession", set())
        self.dict["name"] = prefix
        self.dict["accession"].add(prefix)
        self.dict["organism"] = self.organism


    def add_genome_parameter(self, prefix, value):
        cat = prefix.split(".", 3)[1]
        self.dict.setdefault(cat, set())
        self.dict[cat].add(value)


    def error_logging(self):
        self.error += 1
        f = open(os.path.join(self.currdir, "outputs/errors.txt"), "a")
        f.write(self.dict["name"] + "\n" +
                traceback.format_exc() + "\n" +
                "================================" + "\n")
        print "Error %d occurred." % self.error


    def create_pending_genome(self):
        g = Graph()
        n = self.dict["name"]
        kwargs = {}

        if check_NamedIndividual(n):
            print n + " already in Blazegraph."
        else:
            kwargs.update(self.get_ncbi_ids(n))

            for key, value in self.dict.iteritems():
                if key in self.genome_params:
                    param = self.genome_params.get(key)
                    kwargs.update({param:value})

                elif key == "serotype":
                    kwargs.update(self.get_serotypes(value))

                else:
                    kwargs.update({key:value})

            PendingGenome(g, **kwargs).rdf()
            upload_data(generate_output(g))

    def get_ncbi_ids(self, n):
        bioproject = set()
        biosample = set()

        nuccore_id = return_esearch_uid("nuccore", n)

        bioproject = bioproject | return_elink_uid("nuccore","bioproject",nuccore_id)
        biosample = biosample | return_elink_uid("nuccore", "biosample", nuccore_id)

        return {"bioproject": bioproject, "biosample": biosample}


    def get_serotypes(self, serotypes):
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

        return {"Otype":Otype, "Htype":Htype}

