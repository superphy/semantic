#!/usr/bin/env python

"""Routines for normalizing metadata inputs

Functions take input metadata values and transform them
to standardized values

All functions have the same prototype

Args:
	v(str): Input values

Returns:
	tuple: (
		bool: True if routine made change,
		str: Updated or original value string
		)


"""

import re

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"



def basic_formatting(v):
	"""Formating applied to every value

	Removes trailing and leading whitespace, brackets,
	and transforms to lowercase

	"""
	
	v = re.sub(r'^\s', '', v) # leading whitespace
	v = re.sub(r'\s$', '', v) # trailing whitespace
	v = re.sub(r'\[\]\(\)', '', v) # brackets

	return (True, v.lower());
}