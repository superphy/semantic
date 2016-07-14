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

