import rdflib
import os
import requests
import subprocess
from Bio import SeqIO
import string

from SuperPhy.models.upload.classes.namespaces import n, owl, rdf, xml, xsd, rdfs, gfvo, faldo

class Genomes(object):
    """
    .
    """
    def __init__(self):
        """
        .
        """
        self.data = rdflib.Graph()


    def upload(self):
        """
        post the graph to the triplestore database
        """
        #todo: change this to a call to models/sparql/endpoint
        request = requests.post(
            os.getenv(
                'SUPERPHY_RDF_URL',
                "http://localhost:9000/blazegraph/namespace/superphy/sparql"
            ),
            data=self.data.serialize(format="turtle"),
            headers={'Content-Type':'application/x-turtle'}
        )
        return request.content


    def add_sequence(self, filepath):
        """
        """
        bp = 0
        contigs = 0
        accession = None
        sequences = []

        for record in SeqIO.parse(filepath, "fasta"): #Contigs
            #if self.accession is not None:
            contig = record.id.split('.').pop(0)
            if len(contig) == 12:
                all_ = string.maketrans('', '')
                nochars = all_.translate(all_, string.ascii_letters)
                accession = contig.translate(all_, nochars)
            else:
                accession = contig
            contigs += 1
            bp += len(record.seq)
            sequences.append((">" + accession, record.seq))

        name = accession + "_seq"
        checksum = subprocess.Popen("md5sum {}".format(filepath),\
            shell=True, stdout=subprocess.PIPE).stdout.read(32)
        #is_from = is_from

        self.data.add((
            n[name], rdf.type, owl.NamedIndividual
        ))
        self.data.add((
            n[name], rdf.type, n.Sequence
        ))

        for sequence in sequences:
            self.data.add((
                n[name], n.has_value, rdflib.Literal(str(sequence), datatype=rdflib.XSD.string)
            ))
        self.data.add((
            n[name], n.has_base_pair, rdflib.Literal(str(bp), datatype=rdflib.XSD.string)
        ))
        self.data.add((
            n[name], n.has_contigs, rdflib.Literal(str(contigs), datatype=rdflib.XSD.string)
        ))
        self.data.add((
            n[name], n.has_checksum, rdflib.Literal(str(checksum), datatype=rdflib.XSD.string)
        ))
        """
        #I don't know how to do this properly. - Bryce
        # These literals have historically been attached to this tag:
        # "CORE", 'PLASMID', or "WGS"
        self.data.add((
            n[name], n.is_from, rdflib.Literal(str("WGS"), datatype=rdflib.XSD.string)
        ))
        """
        self.data.add((
            n[accession], n.has_sequence, n[name]
        ))
        self.data.add((
            n[name], n.is_sequence_of, n[accession]
        ))

    def add_metadata(self, data):
        """
        This method is to add data fields to the rdflib Graph object.

        Inputs: A simple key-value dictionary containing at least an
            "Accession" key.
        Changes: The graph object will have triples added to it.
        Outputs: None

        Note: Don't do any validation here. You should be able to make other
        methods that query the rdflib object for appropriatly created metadata
        objects. Do validation between calling this method and uploading the
        genome to an external database. (Or integrate validation with upload.)

        Required Keys:
            Accession

        Optional Keys:
            Serotype_H
            Serotype_O
            Bioproject_id
            Biosample_id
            Isolation_date
            Isolation_location
            Isolation_host
            Isolation_source
            Strain
            Syndrome
            User

        """

        metadata = {}

        for key in data:
            metadata[key.lower()] = data[key]

        #Get the accession # as an unique id.
        name = metadata["accession"]

        #Create a unique object in with name.
        self.data.add((
            n[name], rdf.type, owl.NamedIndividual
        ))

        #Associate the new object with it being a genome
        self.data.add((
            n[name], rdf.type, gfvo.Genome
        ))

        #Add the literal accession value to the object
        self.data.add((
            n[name], n.has_accession, rdflib.Literal(name, datatype=rdflib.XSD.string)
        ))

        #Set the organism of the genome
        self.data.add((
            n[name], n.is_genome_of, n["ecoli"]
        ))
        #Set the object as a genome of E coli.
        #This would not be nessesary if the reasoner is working. This is here
        #   still because it was still in the code I am refactoring.
        self.data.add((n["ecoli"], n.has_genome, n[name]))

        #Set the H serotype to value or 'None'
        try:
            self.data.add((n[name], n.has_Htype, n["H" + str(metadata["serotype_h"])]))
            self.data.add((n["H" + str(metadata["serotype_h"])], n.is_Htype_of, n[name]))
        except KeyError:
            self.data.add((n[name], n.has_Htype, n["H" + str(None)]))
            self.data.add((n["H" + str(None)], n.is_Htype_of, n[name]))

        #Set the O serotype to value or 'None'
        try:
            self.data.add((n[name], n.has_Otype, n["O" + str(metadata["serotype_o"])]))
            self.data.add((n["O" + str(metadata["serotype_o"])], n.is_Otype_of, n[name]))
        except KeyError:
            self.data.add((n[name], n.has_Otype, n["O" + str(None)]))
            self.data.add((n["O" + str(None)], n.is_Otype_of, n[name]))

        #Optional
        #Set the bioproject id metadata
        if "bioproject_id" in metadata:
            self.data.add((
                n[name], n.has_bioproject,\
                rdflib.Literal(metadata["bioproject_id"], datatype=rdflib.XSD.string)
            ))

        #Optional
        #Set the biosample id metadata
        if "biosample_id" in metadata:
            self.data.add((
                n[name], n.has_biosample,\
                rdflib.Literal(metadata["biosample_id"], datatype=rdflib.XSD.string)
            ))

        #Optional
        #Set the isolation date meta-data
        if "isolation_date" in metadata:
            self.data.add((
                n[name], n.has_isolation_date,\
                rdflib.Literal(metadata["isolation_date"], datatype=rdflib.XSD.string)
            ))

        #Optional
        #Set the isolation location metadata
        if "isolation_location" in metadata:
            self.data.add((
                n[name], n.has_geographic_location,\
                rdflib.Literal(metadata["isolation_location"], datatype=rdflib.XSD.string)
            ))
        #Optional
        if "isolation_host" in metadata:
            self.data.add((
                n[name], n.has_isolation_attribute, n[metadata["isolation_host"]]
            ))
            self.data.add((
                n[metadata["isolation_host"]], n.is_isolation_attribute_of, n[name]
            ))

        #Optional
        if "isolation_source" in metadata:
            self.data.add((
                n[name], n.has_isolation_attribute, n[metadata["isolation_source"]]
            ))
            self.data.add((
                n[metadata["isolation_source"]], n.is_isolation_attribute_of, n[name]
            ))

        #Optional
        if "strain" in metadata:
            self.data.add((
                n[name], n.has_strain,\
                rdflib.Literal(metadata["strain"], datatype=rdflib.XSD.string)
            ))

        #Optional
        if "syndrome" in metadata:
            self.data.add((n[name], n.has_isolation_attribute, n[metadata["syndrome"]]))
            self.data.add((n[metadata["syndrome"]], n.is_isolation_attribute_of, n[name]))

        if "user" in metadata:
            self.data.add((
                n[name], n.is_owned_by, n[metadata["user"]]
            ))
            self.data.add((n[metadata["user"]], n.owns, n[name]))
        else:
            self.data.add((
                n[name], n.is_owned_by, n["public"]
            ))
            self.data.add((
                n["public"], n.owns, n[name]
            ))

