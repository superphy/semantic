__author__ = 'Stephen Kan'

from ontology_uploader import upload_all_ontologies
from blazegraph_setup import generate_all

generate_all()
upload_all_ontologies()