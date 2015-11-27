#!/usr/bin/env python

"""Classes for mapping input metadata terms to superphy metadata terms


Classes store metadata and verify correctness and possible conflicts in metadata

Example:
    Usage:

        $ python example_google.py


"""

from pprint import pprint
import json
import logging
from superphy.uploader.external_mapping import cleanup 
from superphy.uploader.external_mapping import validate


__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"


class DecisionTree(object):
	"""Manages instructions for mapping input metadata to superphy ontology



	"""

	def __init__(self, decision_tree_json_file, logger=None):
        """Constructor

        Loads JSON decision tree from file and validates

        Args:
            decision_tree_json_file(str): filepath to decision tree JSON file

        """

        self.logger = logger or logging.getLogger(__name__)

        # Initialize standard parsing routines
        self._default_cleanup_routines = [
            'basic_formatting'
        ]

        # Read json file
 		with open(decision_tree_json_file) as data_file:    
    		self._decision_tree = json.load(data_file)

    	# Validate json file



    def _validate(self):
    	"""Checks decision tree file has correct format and
    	uses defined methods

    	Args:
    		None

    	Returns:
    		bool: True if format is valid

    	"""

    	for term, tree  in self._decision_tree.iteritems():

    		recognized_fields = ['cleanup_routines', 'validation_routines']
    		
    		if tree['keep']:
    			# This term will be parsed


    		else:
    			# This term will be discarded
    			# The block should be empty, if its not warn user

    			for f in recognized_fields:
    				if f in tree:
    					self.logger.warn("Invalid decision tree format: discarded attribute %s should be empty."%term)
    					return False
