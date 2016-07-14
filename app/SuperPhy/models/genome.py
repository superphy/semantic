import rdflib
import os
import requests

from SuperPhy.models.upload.classes.namespaces import n, owl, rdf, xml, xsd, rdfs, gfvo, faldo

class Genome(object):
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

        #Add literal values to the graph
        has_values = {
            "Accession": n.has_accession,
            "Bioproject_Id": n.has_bioproject,
            "Biosample_Id": n.has_biosample,
            "Serotype_H": n.has_Htype,
            "Serotype_O": n.has_Otype,
            "Strain": n.has_strain,
            "Isolation_date": n.has_isolation_date,
            "Isolation_location": n.has_geographic_location,
        }

        for key, predicate in has_values.items():
            if key in metadata:
                self.data.add((
                    n[metadata[key]], predicate, rdflib.Literal(metadata[key], datatype=rdflib.XSD.string)
                ))

        #Set mandatory links
        self.data.add()

genome = Genome()
genome.add_metadata({
    "Accession":"AAJW00000000",
    "Bioproject_Id":"15578",
    "Biosample_Id":"2435896",
    "Serotype_H":"9",
    "Serotype_O":"111",
    "Strain":"E110019"
})
print genome.data.serialize(format="turtle")