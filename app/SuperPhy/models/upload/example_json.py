#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from SuperPhy.models.sparql.genomes import *
from SuperPhy.models.upload._sparql import *

"""
Example using the _sparql functions for returning data
"""

results = find_genome("AYQH01000499")
print results