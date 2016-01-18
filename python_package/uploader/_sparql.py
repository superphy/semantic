#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This module wraps often-used queries to the Blazegraph SPARQL endpoint.
"""

from SPARQLWrapper import JSON, SPARQLWrapper
from superphy import endpoint
from superphy.endpoint import query as _sparql_query
from superphy.endpoint import update as _sparql_update

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"


def find_from_host(host):
    """
    Finds the correct isolation_from_host instance in Blazegraph given a host descriptor, or returns none if nothing
    is found

    Args:
        host: a term used to identify the host (common or scientifi name, generally)

    Returns: the SPARQL URI for the associated isolation_from_host object or None

    """
    results = _sparql_query(
        'PREFIX : <https://github.com/superphy#>\n'
        'SELECT ?p WHERE {?s ?o "%s"^^xsd:string . ?s :is_object_of ?p . ?p rdf:type :isolation_from_host}' % host
    )

    return results["results"]["bindings"][0]["p"]["value"].split("#", 1)[1]


def find_syndrome(syndrome):
    """
    Finds the correct isolation_syndrome instance in Blazegraph given a term, or returns none if nothing
    is found

    Args:
        syndrome: a term used to identify the isolation_syndrome

    Returns: the SPARQL URI for the associated isolation_syndrome or None

    """
    results = _sparql_query(
        'PREFIX : <https://github.com/superphy#>\n'
        'SELECT ?s WHERE {?s ?o "%s"^^xsd:string . ?s rdf:type :isolation_syndrome . }' % syndrome
    )

    return results["results"]["bindings"][0]["s"]["value"].split("#", 1)[1]

def find_source(source):
    """
    Finds the correct isolation_from_source instance in Blazegraph given a term, or returns none if nothing is found

    Args:
        source: a term used to identify the isolation_from_source

    Returns: the SPARQL URI for the associated isolation_from_source or None

    """
    results = _sparql_query(
        'PREFIX : <https://github.com/superphy#>\n'
        'SELECT ?s WHERE {?s ?o "%s"^^xsd:string . ?s rdf:type :isolation_from_source}' % source
    )

    return results["results"]["bindings"][0]["s"]["value"].split("#", 1)[1]


def check_NamedIndividual(name):
    """
    Checks to see if a given SPARQL URI is an instance of any RDF class encoded into the database

    Args:
        name: the SPARQL URI of the instance (must be from the superphy ontology)

    Returns: a boolean indicating if the instance exists or not in the database

    """
    results = _sparql_query(
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX owl: <http://www.w3.org/2002/07/owl#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'ASK { :%s rdf:type owl:NamedIndividual .}' % name
    )

    return results["boolean"]


def find_missing_sequences():
    """
    Finds Genome instances in Blazegraph that are missing a sequence and hasn't failed sequence validation

    Returns:  list of SPARQL URIs for Genome instances

    """
    results = _sparql_query(
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'SELECT ?s ?acc WHERE { ?s rdf:type gfvo:Genome . ?s :has_accession ?acc . '
        'MINUS { ?s :has_valid_sequence ?o }}'
    )

    return ((result["s"]["value"].rsplit("#", 1)[1], result["acc"]["value"])
            for result in results["results"]["bindings"])


def find_duplicate_biosamples():
    """
    Checks to see if a BioSample id is unique or not; if it is not, identify all Genomes that refer to it

    Returns: a list of tuples composed of a BioSample id and a list of SPARQL URIs for Genomes

    """
    results = _sparql_query(
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'SELECT ?BioSample (GROUP_CONCAT( ?Genome ; SEPARATOR = "#") AS ?Genomes) (COUNT (?Genome) AS ?Elements)\n'
        'WHERE { ?Genome rdf:type gfvo:Genome . ?Genome :has_biosample ?BioSample . '
        'MINUS { ?Genome :has_sequence ?Sequence . ?Sequence :is_from "WGS"^^xsd:string .}}\n'
        'GROUP BY ?BioSample HAVING ( ?Elements > 1)'
    )

    return ((result["BioSample"]["value"], result["Genomes"]["value"].split("#", )[1::2])
            for result in results["results"]["bindings"])


def find_core_genome(biosample):
    """
    Finds all Genomes with the specified BioSample id that are core genomes (labelled with CORE)

    Args:
        biosample: BioSample id of interest

    Returns: a list of SPARQL URIs of Genomes that match the BioSample and are core genomes

    """
    results = _sparql_query(
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'SELECT ?Genome \n'
        'WHERE { ?Genome rdf:type gfvo:Genome . ?Genome :has_biosample "%s"^^xsd:string. '
        '?Genome :has_sequence ?Sequence . ?Sequence :is_from "CORE"^^xsd:string .}\n' % biosample
    )

    return [result["Genome"]["value"].split("#", 1)[1] for result in results["results"]["bindings"]]


def find_category(category):
    """
    Finds the corresponding category instance in Blazegraph. Returns None if nothing is found.

    Args:
        category: a term used to identify the category

    Returns: the SPARQL URI for the associated category
    """
    results = _sparql_query(
        'PREFIX : <https://github.com/superphy#>\n'
        'SELECT ?s WHERE {?s ?o "%s"^^xsd:string . ?s rdf:type :sub_category}' % category 
    )

    return results["results"]["bindings"]


def find_sub_category(sub_category):
    """
    Finds the corresponding sub_category instance in Blazegraph. Returns None if nothing is found.

    Args:
        sub_category: a term used to identify the sub_category

    Returns: the SPARQL URI for the associated sub_category
    """
    results = _sparql_query(
        'PREFIX : <https://github.com/superphy#>\n'
        'SELECT ?s WHERE {?s ?o "%s"^^xsd:string . ?s rdf:type :sub_category}' % sub_category 
    )

    return results["results"]["bindings"]


def find_genome(accession):
    """
    Finds the genome instance in Blazegraph. Returns None if nothing is found.

    Args:
        genome: genome accession number

    Returns: the SPARQL URI for the associated genome instance. Returns None if nothing found.
    """
    results = _sparql_query(
        'PREFIX : <https://github.com/superphy#>\n'
        'SELECT ?s WHERE {?s :has_accession "%s"^^xsd:string . }' % accession 
    )

    return results["results"]["bindings"][0]["s"]["value"]


def delete_instance(name):
    """
    Deletes an instance with a given SPARQL URI on the database by removing all triples with it (assumption:
    not a predicate, but then predicates aren't instances)

    Args:
        name: the SPARQL URI of the instance you want to delete

    Prints out the response from the server regarding the SPARQL Update query

    """
    print _sparql_update(
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'DELETE { :%s ?property ?object . ?subject ?property :%s . }\n'
        'WHERE { { :%s ?property ?object } UNION { ?subject ?property :%s } }' % (name, name, name, name)
    )

def insert_accession_sequence(core, plasmid, plasmid_seq):
    """
    Given the SPARQL URIs for the core genome, the plasmid genome, and plasmid sequence, adds the plasmid sequence
    under the core genome and removes the connection to the plasmid genome.

    This is a cleanup routine, as core and plasmid genomes would share the same metadata, and it would not make sense
    to have both a plasmid and core genome instance containing the same metadata.

    Args:
        core (str): SPARQL URI based on the superphy ontology for the core genome instance
        plasmid (str): SPARQL URI based on the superphy ontology for the plasmid genome instance
        plasmid_seq (str): SPARQL URI based on the superphy ontology for the plasmid sequence instance
\
    Prints out the response from the server regarding the SPARQL Update query

    """
    print _sparql_update(
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'INSERT DATA { :%s :has_accession "%s"^^xsd:string. '
        ':%s :has_sequence :%s . :%s :is_sequence_of :%s . }\n'
        % (core, plasmid, core, plasmid_seq, plasmid_seq, core)
    )


def check_blank_nodes():
    """
    Checks to see if there are any blank nodes present on the database

    Returns: a boolean indicating if blank nodes exists or not in the database

    """
    results = _sparql_query(
        'ASK {?x ?y ?z . FILTER ( isBlank(?x) || isBlank(?z) )}'
    )

    return results["boolean"]


def delete_blank_nodes():
    """
    Deletes all blank nodes on Blazegraph.

    This should only be ran when setting up the database, as blank nodes increase exponentially as triples get
    added and will cause a memory overflow error. Many of these blank nodes are from the interactions of the ontologies
    that superphy is built upon, but they do not significantly contribute to querying.

    Prints out the response from the server regarding the SPARQL Update query

    """
    print _sparql_update(
        'DELETE { ?x ?y ?z }'
        'WHERE { ?x ?y ?z . FILTER ( isBlank(?x) || isBlank(?z) ) }'
    )


def check_checksum(checksum):
    """
    Checks if a particular checksum exists in the database.

    As checksums are supposed to be unique to the sequence, if any are found in the database, the chances of
    there being a duplicate sequence is high.

    Args:
        checksum (str): the hash for a sequence

    Returns: a boolean indicating if the hash was found in the database

    """
    results = _sparql_query(
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'ASK { ?Sequence :has_checksum "%s"^^xsd:string}' % checksum
    )

    return results["boolean"]


'''
def _sparql_query(query):
    sparql = SPARQLWrapper("http://localhost:9999/blazegraph/namespace/superphy/sparql")
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results

def _sparql_update(query):
    sparql = SPARQLWrapper("http://localhost:9999/blazegraph/namespace/superphy/sparql")
    sparql.method = 'POST'
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results
'''
