#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""

"""

import os
import subprocess

from SuperPhy.models.upload._utils import generate_path

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
