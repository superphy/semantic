__author__ = 'Stephen Kan'

from ontology_uploader import upload_all_ontologies, create_namespace
from blazegraph_setup import generate_all
import gc

create_namespace()
generate_all()
upload_all_ontologies()
gc.collect()