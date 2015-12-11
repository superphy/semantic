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


