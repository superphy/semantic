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
        .
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

        Inputs: A simple key-value dictionary containing values for
        Changes: The graph object will have triples added to it.
        Outputs: None
        """
        #'accession', 'bioproject', 'biosample', 'htype', 'isolation_date', 'isolation_host', 'isolation_location', 'isolation_source', 'organism', 'otype', 'strain', 'syndrome'

        if "Accession" in metadata:
            self.data.add((
                n[metadata["Accession"]], n.has_accession,\
                rdflib.Literal(metadata["Accession"], datatype=rdflib.XSD.string)
            ))

        if "Bioproject_Id" in metadata:
            self.data.add((
                n[metadata["Bioproject_Id"]], n.has_bioproject,\
                rdflib.Literal(metadata["Bioproject_Id"], datatype=rdflib.XSD.string)
            ))
            
        if "Biosample_Id" in metadata:
            self.data.add((
                n[metadata["Biosample_Id"]], n.has_biosample,\
                rdflib.Literal(metadata["Biosample_Id"], datatype=rdflib.XSD.string)
            ))

        if "Serotype_H" in metadata:
            self.data.add((
                n[metadata["Serotype_H"]], n.has_Htype,\
                rdflib.Literal(metadata["Serotype_H"], datatype=rdflib.XSD.string)
            ))

        if "Serotype_O" in metadata:
            self.data.add((
                n[metadata["Serotype_O"]], n.has_Otype,\
                rdflib.Literal(metadata["Serotype_O"], datatype=rdflib.XSD.string)
            ))

        if "Strain" in metadata:
            self.data.add((
                n[metadata["Strain"]], n.has_strain,\
                rdflib.Literal(metadata["Strain"], datatype=rdflib.XSD.string)
            ))

        if "Isolation_date" in metadata:
            self.data.add((
                n[metadata["Isolation_date"]], n.has_isolation_date,\
                rdflib.Literal(metadata["Isolation_date"], datatype=rdflib.XSD.string)
            ))

        if "Isolation_location" in metadata:
            self.data.add((
                n[metadata["Isolation_location"]], n.has_geographic_location,\
                rdflib.Literal(metadata["Isolation_location"], datatype=rdflib.XSD.string)
            ))

genomes = Genomes()
genomes.add_metadata({
    "Accession":"AAJW00000000",
    "Bioproject_Id":"15578",
    "Biosample_Id":"2435896",
    "Serotype_H":"9",
    "Serotype_O":"111",
    "Strain":"E110019"
})
print genomes.data.serialize(format="turtle")