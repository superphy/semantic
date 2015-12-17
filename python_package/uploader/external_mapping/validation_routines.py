#!/usr/bin/env python

"""Container for all the validation_routine functions

Validation object contains functions to convert foreign attribute-value pairs to Superphy meta-data terms and standard
values.

Validation routine methods have specification:
    Args: 
      v(str): Input string 
      mapper(Mapper): a reference to the mapper object


    Returns:
      False when no match was found
        -or-
      String 'skip' to ignore values like 'NA' or 'missing'
        -or-
      List of tuples containing:
        [(superphy_term_uri, superphy_value_uri),...]
      

"""

import re
from datetime import datetime
from pprint import pprint

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"


##################
# Utility methods
##################

def _exact_match(v, patterns):

    combined = "(" + ")|(".join(patterns) + ")"
    
    if re.match(combined, v):
        return True
    else:
        return False


##################
# Matching methods
##################


def non_value(v, mapper=None):
    """Non-values that should be ignored like 'missing'

    """

    inputs = [
        'not determined',
        'not applicable',
        'not collected',
        'missing',
        'N/A',
        'NA',
        'None',
        'Unknown',
        '0'
    ]
   
    if _exact_match(v, inputs):
        return 'skip'

    else:
        return False


# Hosts
def hosts(v, mapper):
    """Checks for match against known host values

    Note: the cleanup routines should handle any
    synonyms
    
    """

    if _exact_match(v, mapper.ontology('host').uri_list()):
        return [('host', v)]
    else:
        return False


def sources(v, mapper):
    """Checks for match against known source values

    Note: the cleanup routines should handle any
    synonyms
    
    """

    if _exact_match(v, mapper.ontology('source').uri_list()):
        return [('source', v)]
    else:
        return False


def syndromes(v, mapper):
    """Checks for match against known syndrome values

    Note: the cleanup routines should handle any
    synonyms
    
    """

    if _exact_match(v, mapper.ontology('syndrome').uri_list()):
        return [('syndrome', v)]
    else:
        return False


def multi_syndromes(v, mapper):
    """Splits disease descriptions into individual syndrome values
    
    """
  
    breakdowns = {
        'uti induced bacteremia': [
            'uti',
            'bacteremia'
        ],
        'hc, hus': [
            'hc',
            'hus'
        ]
    }
    
    if _exact_match(v, breakdowns.keys ):
        return [ ('syndrome', i) for i in breakdowns[v] ]
    else:
        return False



def serotype(v, mapper):
    """Converts serotype into two-part meta-data values
    
    """

    sero = []
    
    # O-antigen
    oresult = re.match(r"^o([\d+|nt|na])\:", v)

    if oresult:
        match = oresult.group(1)
        if match != 'na':
            sero.append(('o_serotype', match))

    # H-antigen
    hresult = re.match(r"\:h([\d+|nm|na])$", v)

    if hresult:
        match = hresult.group(1)

        if match != 'na':
            sero.append(('h_serotype', match))
    

def collection_date(v, mapper):
    """Store valid date in past

    """

    try:
        dt = datetime.strptime(v, "%Y-%m-%d")

        if dt >= datetime.today():
            return False

        else:
            return ('isolation_date', dt.strftime("%Y-%m-%d"))

    except ValueError:
        return False


def location(v, mapper):
    """Use google maps V3 geocoder
    to validate locations.  Add new valid locations
    to superphy cache

    """

    pass

"""serotype => {
        cvterm => 'serotype',
        priority => 0
    },
    serovar => {
        cvterm => 'serotype',
        priority => 0
    },
    strain => {
        cvterm => 'strain',
        priority => 0
    },
    sub_strain => {
        cvterm => 'strain',
        priority => 1
    },
    culture_collection => {
        cvterm => 'strain',
        priority => 1
    },
    sub_species => {
        cvterm => 'strain',
        priority => 1
    },
    isolate => {
        cvterm => 'strain',
        priority => 1
    },
    collection_date => {
        cvterm => 'isolation_date',
        priority => 0
    },
    pubmed => {
        cvterm => 'pmid',
        priority => 0,
    },
    direct_submission => {
        cvterm => 'owner',
        priority => 0,
    },
    comment => {
        cvterm => 'comment',
        priority => 0
    },
    note => {
        cvterm => 'comment',
        priority => 1
    },
    keyword => {
        cvterm => 'keywords',
        priority => 0
    },
    mol_type => {
        cvterm => 'mol_type',
        priority => 0, 
        static => 1
    },
    finished => {
        cvterm => 'finished',
        priority => 0,
        static => 1
    },
    description => {
        cvterm => 'description',
        priority => 0,
    },
    secondary_dbxref => {
        cvterm => 'secondary_dbxref',
        priority => 0,
    },
    primary_dbxref => {
        cvterm => 'primary_dbxref',
        priority => 0,
        static => 1
    },
    country => {
        cvterm => 'isolation_location',
        priority => 0,
    },
    host => {
        cvterm => 'isolation_host',
        priority => 0,
    },
    isolation_source => {
        cvterm => 'isolation_source',
        priority => 0,
    },
    syndrome => {
        cvterm => 'syndrome',
        priority => 0,
    },
    name => {
        cvterm => 'name',
        priority => 0,
    },
    uniquename => {
        cvterm => 'uniquename',
        priority => 0,
    }

    location (lat,lon) or formatted address as uri
    properties (search_query, raw json, lat,lon, formatted address)
   
"""
