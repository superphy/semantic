#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""

"""

import os
import subprocess

from SuperPhy.models.upload._utils import generate_path

__author__ = 'Stephen Kan'
__copyright__ = """
    © Copyright Government of Canada 2012-2015. Funded by the Government of
    Canada Genomics Research and Development Initiative"""
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'

"""
    The functions in this are commented out, because it doesn't work to be
    copying deleting Gigs of data if you are testing locally. Another solution
    would be to delete the triples after creating them, or restarting the
    client from a different position, forceing the jar file to make a new file
    in a different spot.
"""

class BlazegraphIntegration(object):
    @classmethod
    def setUpClass(cls):
        pass
        
    @classmethod
    def tearDownClass(cls):
        pass
        #top_dir = generate_path("../../../../../../..")
        #os.chdir(top_dir)
        #src = os.path.join(
        #    os.getcwd(),
        #    "superphy/database/blazegraph/testing/bigdata.jnl"
        #)
        #subprocess.call("rm -f %s" % src, shell=True)
