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
from sparql import check_NamedIndividual

reload(sys)
sys.setdefaultencoding("utf-8")


class MinerDataUploader(object):
    genome_params = {"isolation_date":"date", "isolation_location":"location", "isolation_host":"host",
                     "isolation_source":"source"}
    currdir = os.path.dirname(inspect.getfile(inspect.currentframe()))

    def __init__(self, filename, organism):
        self.progress = 0
        self.error = 0
        self.dict = {}
        self.filename = filename
        self.organism = organism

        self.iterate()


    def load_JSON(self):
        path = os.path.join(self.currdir, self.filename)
        fd = open(path, "r")
        self.parser = ijson.parse(fd)


    def iterate(self):
        self.load_JSON()
        for prefix, event, value in self.parser:
            if ("." not in prefix) and (prefix is not "") and (event == "start_map"):
                self.add_to_graph()
                self.start_new_genome(prefix)

            if prefix.endswith(".displayname"):
                self.add_genome_parameter(prefix, value)

        self.add_to_graph()
        print "%d genomes parsed, %d errors occurred." %(self.progress, self.error)


    def start_new_genome(self, prefix):
        self.dict.setdefault("accession", set())
        self.dict["name"] = prefix
        self.dict["accession"].add(prefix)
        self.dict["organism"] = self.organism


    def add_genome_parameter(self, prefix, value):
        cat = prefix.split(".", 3)[1]
        self.dict.setdefault(cat, set())
        self.dict[cat].add(value)


    def add_to_graph(self):
        if self.dict:
            try:
                self.progress += 1
                print str(self.progress) + ": downloading files"
                self.create_pending_genome()
            except Exception as e:
                self.error_logging()

            self.dict.clear()


    def error_logging(self):
        f = open(os.path.join(self.currdir, "outputs/errors.txt"), "a")
        f.write(self.dict["name"] + "\n")
        f.write(traceback.format_exc() + "\n" + "==============================================" + "\n")
        self.error += 1
        print "Error %d occurred." % self.error


    def create_pending_genome(self):
        g = Graph()
        n = self.dict["name"]
        kwargs = {}

        if check_NamedIndividual(n):
            print n + " already in Blazegraph."
        else:
            self.add_ids(kwargs, n)

            for key, value in self.dict.iteritems():
                if key in self.genome_params:
                    param = self.genome_params.get(key)
                    kwargs.update({param:value})

                elif key == "serotype":
                    (Otype, Htype) = self.return_serotypes(value)
                    kwargs.update({"Otype":Otype, "Htype":Htype})

                else:
                    pass
                    kwargs.update({key:value})

            PendingGenome(g, **kwargs).rdf()
            ontology_uploader.upload_data(generate_output(g))


    def add_ids(self, kwargs, n):
        nuccore_id = return_esearch_uid("nuccore", n)
        (biosample, bioproject) = self.return_bio_ids(nuccore_id)
        kwargs.update({"bioproject": bioproject, "biosample": biosample})


    def return_serotypes(self, serotypes):
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


    def return_bio_ids(self, nuccore_id):
        bioproject = set()
        biosample = set()

        try:
            bioproject = bioproject | return_elink_uid("nuccore","bioproject",nuccore_id)
            biosample = biosample | self.return_elink_biosample(bioproject, biosample, nuccore_id)

        except IndexError:
            self.return_efetch(bioproject, biosample, nuccore_id)

        return (bioproject, biosample)


    def return_efetch(self, bioproject, biosample, nuccore_id):
        for record in return_nuccore_efetch(nuccore_id):
            for xref in record["GBSeq_xrefs"]:
                if xref["GBXref_dbname"] == "BioProject":
                    bioproject.add(only_digits(xref["GBXref_id"]).lstrip("0"))
                elif xref["GBXref_dbname"] == "BioSample":
                    biosample.add(only_digits(xref["GBXref_id"]).lstrip("0"))


    def return_elink_biosample(self, bioproject, biosample, nuccore_id):
        item = set()
        try:
            item = biosample | return_elink_uid("nuccore", "biosample", nuccore_id)

        except IndexError:
            for id in bioproject:
                item = biosample | return_elink_uid("bioproject", "biosample", id)
        return item