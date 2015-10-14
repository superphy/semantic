__author__ = 'ubiquitin'

import blazegraph_setup
import ontology_uploader

ontology_uploader.upload_all_ontologies()
ontology_uploader.upload_ontology("outputs/results.ttl")


