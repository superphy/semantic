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
from pprint import pprint

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
	v = str(v) # Convert to byte string
	v = re.sub(r'^\s', '', v) # leading whitespace
	v = re.sub(r'\s$', '', v) # trailing whitespace
	v = re.sub(r'\[\]\(\)', '', v) # brackets

	return (True, v.lower())


def _replacement(v, patterns, repl):
	"""If value matches any item in list, return replacement string

	"""

	combined = "(" + ")|(".join(patterns) + ")"

	clean_v, num = re.subn(combined, repl, v, flags=re.I)

	if num > 0:
		return (True, clean_v)
	else:
		return (False, v)


def fix_hosts(v):
	"""Runs individual host methods until recognized host is detected

	For this function to work properly, there can be no overlapping terms
    between routines

	"""

	routines = [
		'fix_human',
		'fix_cow'
	]

	success = False
	for r in routines:
		routine = globals()[r]

		cleaned, clean_v = routine(v)
		if cleaned:
			v = clean_v
			success = True
			break

	return(success, v)


def fix_human(v):
	"""Map human synonyms to human ontology uri

	"""

	synonyms = [
		'Human Homo sapiens',
		'Homo sapiens',
		'human',
		'patient',
		'infant',
		'child',
		'9606'
	]

	return _replacement(v, synonyms, 'hsapiens')


def fix_cow(v):
	"""Map cow synonyms to cow ontology uri

	"""

	synonyms = [
		"Bos taurus",
		"cow",
		"cattle",
		"calf",
		"bovine",
	]

	return _replacement(v, synonyms, 'btaurus');


def fix_sources(v):
	"""Runs individual source methods until recognized source is detected

	For this function to work properly, there can be no overlapping terms
    between routines

	"""

	routines = [
		'fix_poop',
		'fix_intestine',
		'fix_ut'
	]

	success = False
	for r in routines:
		routine = globals()[r]

		cleaned, clean_v = routine(v)
		if cleaned:
			v = clean_v
			success = True
			break

	return(success, v)


def fix_poop(v):
	"""Map feces synonyms to feces ontology uri

	"""

	synonyms = [
		'stool sample',
		'fecal sample',
		'feces envo:00002003',
		'animal - feces',
		'stool',
		'fecal',
		'feces'
	]

	return _replacement(v, synonyms, 'feces')


def fix_intestine(v):
	"""Map intestine synonyms to intestine ontology uri

	"""

	synonyms = [
		'stool sample',
		'fecal sample',
		'feces envo:00002003',
		'animal - feces',
		'stool',
		'fecal',
		'feces'
	]

	return _replacement(v, synonyms, 'feces')
