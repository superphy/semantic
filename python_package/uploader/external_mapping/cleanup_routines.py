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
from dateutil.parser import parse
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
		"gastrointestinal_tract", 
		"intestinal mucosa tissue",
		"gastrointestinal"
	]

	return _replacement(v, synonyms, 'intestine')


# Convert all urinary tract synonyms
def fix_ut(v):
	"""Map intestine synonyms to intestine ontology uri

	"""

	synonyms = [
		"urinary tract",
		"urogenital_tract",
		"genitourinary"
	]
	
	return _replacement(v, synonyms, 'urogenital')


def fix_water(v):
	"""Map water source synonyms to environment:water uri

	"""
	
	synonyms = [
		"terrain - watershed",
		"water - canal",
		"water - river",
		"water - stream",
		"water - intake",
		"water - waste water",
		"agricultural - irrigation ditch",
		"environmental water study",
		"water study"
		"water"
	]
	
	return _replacement(v, synonyms, 'environment: water');


def fix_syndromes(v):
	"""Map several disease synonyms to their DB uri

	"""
	
	success = False

	diseases = {
		"uti": ['urinary tract infection', 'recurrent uti'],
		"hus": ['hemolytic uremic syndrome'],
		"hc": ['hemorrhagic colitis'],
		"septicaemia": ['sepsis'],
		"diarrhea": ['travellers diarhhea']
	}
		
	for d in diseases:
		cleaned, clean_v = _replacement(v, diseases[d], d)

		if cleaned: 
			v = clean_v
			success = True
		
	return (success, v);


def fix_serotype(v):
	"""Edit serotype so the follow standard nomenclature

	"""

	# Run serotype through some regex 'cleaners';
	
	# O fixes
	new_v = re.sub(r'\:k\d+', '', v, flags=re.IGNORECASE) # Remove capsule
	new_v = re.sub(r'\s*non-typable', 'nt', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'^ountypeable', 'ont', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'^o\?', 'ont', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'e. coli\s*\b', '', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'sf\?', '', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'^or', 'ont', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'^0', 'o', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'^(\d)', r'o\1', new_v, flags=re.IGNORECASE) # Put o in front of leading number
	new_v = re.sub(r'^(o\d+)[a-z]+\:', r'\1\:', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'^o\d+\s*\,\s*o\d+[a-z]*+\:', r'ont\:', new_v, flags=re.IGNORECASE)
	

	# H fixes
	new_v = re.sub(r'\:h-$', r'\:hnm', new_v, flags=re.IGNORECASE) # Missing H
	new_v = re.sub(r'\:hnm$', r'\:hnm', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'\:ut$', r'\:hnm', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'\:h\s*rough$', r'\:hnm', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'-$', r'\:hnm', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'\:huntypeable$', r'\:hnm', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'\:\?$', r'\:hna', new_v, flags=re.IGNORECASE)
	new_v = re.sub(r'(o\d+)$', r'\1\:hna', new_v, flags=re.IGNORECASE)
	
	# Remove trailing ???
	new_v = re.sub(r'\s*\?+$', r'', new_v, flags=re.IGNORECASE)
	
	if new_v != v:
		return (True, new_v)
	else:
		return (False, v)


def parse_date(v):
	"""Try to parse date and convert to standard format

	"""

	try:
		datetime = parse(v)
		return (True, datetime.strftime('%Y-%m-%d'))
	except ValueError:
		return False
	
