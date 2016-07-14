import rdflib
import os
import requests

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

    def add_metadata(self, metadata):
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

        #Get the accession # as an unique id.
        name = metadata["Accession"]

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
            self.data.add((n[name], n.has_Htype, n["H" + str(metadata["Serotype_H"])]))
            self.data.add((n["H" + str(metadata["Serotype_H"])], n.is_Htype_of, n[name]))
        except KeyError:
            self.data.add((n[name], n.has_Htype, n["H" + str(None)]))
            self.data.add((n["H" + str(None)], n.is_Htype_of, n[name]))

        #Set the O serotype to value or 'None'
        try:
            self.data.add((n[name], n.has_Otype, n["O" + str(metadata["Serotype_O"])]))
            self.data.add((n["O" + str(metadata["Serotype_O"])], n.is_Otype_of, n[name]))
        except KeyError:
            self.data.add((n[name], n.has_Otype, n["O" + str(None)]))
            self.data.add((n["O" + str(None)], n.is_Otype_of, n[name]))

        #Optional
        #Set the bioproject id metadata
        if "Bioproject_id" in metadata:
            self.data.add((
                n[name], n.has_bioproject,\
                rdflib.Literal(metadata["Bioproject_id"], datatype=rdflib.XSD.string)
            ))

        #Optional
        #Set the biosample id metadata
        if "Biosample_id" in metadata:
            self.data.add((
                n[name], n.has_biosample,\
                rdflib.Literal(metadata["Biosample_id"], datatype=rdflib.XSD.string)
            ))

        #Optional
        #Set the isolation date meta-data
        if "Isolation_date" in metadata:
            self.data.add((
                n[name], n.has_isolation_date,\
                rdflib.Literal(metadata["Isolation_date"], datatype=rdflib.XSD.string)
            ))

        #Optional
        #Set the isolation location metadata
        if "Isolation_location" in metadata:
            self.data.add((
                n[name], n.has_geographic_location,\
                rdflib.Literal(metadata["Isolation_location"], datatype=rdflib.XSD.string)
            ))
        #Optional
        if "Isolation_host" in metadata:
            self.data.add((
                n[name], n.has_isolation_attribute, n[metadata["Isolation_host"]]
            ))
            self.data.add((
                n[metadata["Isolation_host"]], n.is_isolation_attribute_of, n[name]
            ))

        #Optional
        if "Isolation_source" in metadata:
            self.data.add((
                n[name], n.has_isolation_attribute, n[metadata["Isolation_source"]]
            ))
            self.data.add((
                n[metadata["Isolation_source"]], n.is_isolation_attribute_of, n[name]
            ))

        #Optional
        if "Strain" in metadata:
            self.data.add((
                n[name], n.has_strain,\
                rdflib.Literal(metadata["Strain"], datatype=rdflib.XSD.string)
            ))

        #Optional
        if "Syndrome" in metadata:
            self.data.add((n[name], n.has_isolation_attribute, n[metadata["Syndrome"]]))
            self.data.add((n[metadata["Syndrome"]], n.is_isolation_attribute_of, n[name]))

        if "User" in metadata:
            self.data.add((
                n[name], n.is_owned_by, n[metadata["User"]]
            ))
            self.data.add((n[metadata["User"]], n.owns, n[name]))
        else:
            self.data.add((
                n[name], n.is_owned_by, n["public"]
            ))
            self.data.add((
                n["public"], n.owns, n[name]
            ))

