#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This module queries NCBI using the Entrez E-utilities programming tools provided by a BioPython wrapper.

As many queries follow the same format and uses the same database, it is expedient to wrap them into functions
requiring few parameters.

References:
    Entrez Programming Utilities Help (http://www.ncbi.nlm.nih.gov/books/NBK25501/)
        - the section "The E-utilities In-Depth: Parameters, Syntax and More" explains how to customize requests
    BioPython documentation (http://biopython.org/DIST/docs/api/Bio.Entrez-module.html)
    BIoPython Tutorial (http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc108)
"""

from Bio import Entrez

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

Entrez.email = "superphy.info@gmail.com"

def return_esearch_uid(db, term):
    """
    Finds the id(s) of entries matching the term from the specified database. This is used to find the gid
    of Nucleotide (db=nuccore) entries from their accession ids (expect an unique match)

    Args:
        db (str): the keyword for the relevant NCBI database
        term (str): the search term used to query the database

    Returns: the first database id found

    """
    handle = Entrez.esearch(db=db, retmax=5, term=term)
    record = Entrez.read(handle)
    return str(record["IdList"][0])

def return_elink_uid(dbfrom, db, id):
    """
    Finds entries in db for an entry in dbfrom that is identified by id

    Args:
        dbfrom (str): the originating database
        db (str): the database with entries of interest
        id (str): the id of an entry in dbform

    Returns: a set of ids for entries in db

    """
    handle = Entrez.elink(dbfrom=dbfrom, db=db, id=id)
    records = Entrez.parse(handle)
    ret = set()

    for record in records:
        ret.add(record["LinkSetDb"][0]["Link"][0]["Id"])

    return ret