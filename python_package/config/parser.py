#!/usr/bin/env python

"""Loads ENV config into python dictionary

Superphy config properties are provided as bash environment variables.  This script 
checks and loads the ENV config settings into a python dictionary. If a 
required property is missing / uninitialized it throws an error.

Example:
        $ from python_package import config.parser.get
        $ config = config.parser.get()

"""

from os import environ

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "CC-SY"
__version__ = "4.0"
__maintainer__ = "Matt"
__email__ = "matthew.whiteside@canada.ca"



env_vars = { 
    'rdf_url': 'SUPERPHY_RDF_URL'
}
"""dictionary: Required config properties and their corresponding ENV variable names

"""

def read(test=None):
    """Reads values from environment variables in 'env_vars' dictionary
    and stores them in a dictionary. 

    Returns:
        dictionary: 
            {
                'env_var_name': value,
            }
        Dictionary keys match key names in 'env_vars' dictionary

    Raises:
        SuperphyConfigError: when variable in 'env_vars' values
        is not defined

    """

    config = {}

    #pprint.pprint(env_vars)

    for prop,var in env_vars.iteritems():
        val = environ.get(var)
        if not val:
            msg = "Superphy environment variable {} is not defined".format(var)
            raise SuperphyConfigError(msg)
        else:
            config[prop] = val

    return config
   

class SuperphyConfigError(Exception):
    """Exception class for missing/uninitialized config.

    """
    pass


