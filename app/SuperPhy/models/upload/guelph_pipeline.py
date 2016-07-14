
import os
import json
import subprocess
import tablib
from Bio import SeqIO
import rdflib

from SuperPhy.models.upload.classes.sequence import Sequence
from SuperPhy.models.upload.blazegraph_upload import BlazegraphUploader
from SuperPhy.models.upload.metadata_upload import MetadataUploader, Metadata, GenomeMetadataUploader, GenomeMetadata
from SuperPhy.models.upload.contig_upload import ContigUploader, ContigsWrapper

class Pipeline(object):
    def __init__(self, fasta_files=list, metadata=str):
        print "FASTA FILES: ", fasta_files, "METADATA: ", metadata
        self.fasta_files = []
        for file_ in fasta_files:
            self.__add_fasta_file(file_)
        self.meta_files = []
        self.__add_csv(metadata)

    def __add_fasta_file(self, file_):
        self.fasta_files.append(file_)

    def __add_csv(self, metadata):
        self.meta_files += self.__split_csv(metadata)

    @classmethod
    def __split_csv(cls, filepath):
        basedir = os.path.realpath(os.path.dirname(__file__)).rsplit("SuperPhy", 1)[0]
        uploads = os.path.join(basedir, 'uploads')

        data = tablib.Dataset()
        data.csv = open(filepath).read()
        print data
        def __validate_headers(data):
            pass
        def __validate_fields(data):
            pass

        __validate_headers(data)
        __validate_fields(data)

        json_files = []
        for genome in data.dict:
            file_ = os.path.join(uploads, genome["Accession"]+".json")
            json_files.append(file_)
            for key, value in genome.items():
                if value == "" or value == []:
                    del genome[key]
            with open(json_files[-1], 'w') as outfile:
                json.dump(genome, outfile, indent=4, sort_keys=True, separators=(',', ':'))
        return json_files

    def __validate_locations(self, files):
        pass

    def process(self):
        output = []

        for file_ in self.fasta_files:
            s = SequenceData()
            s.load(file_)
            #m = GenomeMetadataUploader(file_, "Eschericia Coli")
            #m.upload() #Import, not upload
            s.send_fasta_data()
            output.append(s.show_fasta_data())
        return output



class SequenceData(object):
    def __init__(self):
        self.md5sum = None
        self.base_pairs = 0
        self.contig_numbers = 0
        self.data = []
        self.accession = None
        self.graph = None

    def validate(self):
        """
        Validate that the data is proper.
        """
        pass

    def load(self, filepath):
        """
        Loads the sequence data from a fasta file.
        """
        self.md5sum = subprocess.Popen("md5sum {}".format(filepath),\
            shell=True, stdout=subprocess.PIPE).stdout.read(32)
        self.base_pairs = 0
        self.contig_numbers = 0
        print "FILEPATH: ", filepath
        for record in SeqIO.parse(filepath, "fasta"): #Contigs
            #if self.accession is not None:
            self.accession = ContigsWrapper.genome_name(record.id.split('.').pop(0))
            self.contig_numbers += 1
            self.base_pairs += len(record.seq)
            self.data.append((">" + self.accession, record.seq))

    def send_fasta_data(self):
        """
        Uploads and sends the fasta data to Blazegraph
        """
        self.graph = rdflib.Graph()
        seq = Sequence(self.graph, "{}{}".format(self.accession, "_seq"),\
            self.accession, self.data, self.base_pairs, self.contig_numbers,\
            self.md5sum, "WGS")
        seq.rdf()

        output = self.graph.serialize(format='turtle')

        uploader = BlazegraphUploader
        uploader.upload_data(output)

    def show_fasta_data(self):
        """
        ### This is to return the triples in a human readable format
        ### If you aren't using this, you should remove it.

        """
        data = {}

        # go through all the triples in the graph
        for subject, predicate, object_ in self.graph:
            # initialize property dictionary if we've got a new subject
            if not data.has_key(subject):
                data[subject] = {}
            # initialize object list if we've got a new subject-property combo
            if not data[subject].has_key(predicate):
                data[subject][predicate] = []
            # determine the value dictionary for the object
            value = {'value': unicode(object_)}
            if isinstance(object_, rdflib.URIRef):
                value['type'] = 'uri'
            elif isinstance(object_, rdflib.BNode):
                value['type'] = 'bnode'
            elif isinstance(object_, rdflib.Literal):
                value['type'] = 'literal'
                if object_.language:
                    value['lang'] = object_.language
                if object_.datatype:
                    value['datatype'] = unicode(object_.datatype)
            # add the triple
            data[subject][predicate].append(value)
        return data
