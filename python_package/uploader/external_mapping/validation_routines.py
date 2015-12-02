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

class IV(object):
    """Input validation

    Checks if input is recognized and returns equivalent superphy metatdata term
    and value

    """

    def __init__(self, ontology_dict):
        """


        Args:
            logger (Optional[object]): Pointer to logger op
            param3 (List[str]): Description of `param3`.

        """



##################
# Utility methods
##################

def _exact_match(v, patterns):

    pprint(patterns)

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


