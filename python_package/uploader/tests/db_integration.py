#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""

"""

import os
import subprocess

from superphy.uploader._utils import generate_path

__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'

class BlazegraphIntegration(object):
    @classmethod
    def setUpClass(cls):
        top_dir = generate_path("../../../")
        os.chdir(top_dir)
        src = os.path.join(os.getcwd(),"db/bigdata.jnl")
        dst = os.path.join(os.getcwd(),"db/bigdata.jnl.bk")
        subprocess.call("bash bash/kill_port 9999", shell=True)
        print "Killing existing Blazegraph process"
        subprocess.call("cp %s %s" %(src, dst), shell=True)
        subprocess.call("bash bash/start_blazegraph", shell=True)

    @classmethod
    def tearDownClass(cls):
        top_dir = generate_path("../../../")
        os.chdir(top_dir)
        src = os.path.join(os.getcwd(),"db/bigdata.jnl.bk")
        dst = os.path.join(os.getcwd(),"db/bigdata.jnl")
        subprocess.call("bash bash/kill_port 9999", shell=True)
        print "Killing existing Blazegraph process"
        subprocess.call("cp %s %s" %(src, dst), shell=True)
        subprocess.call("rm -f %s" % src, shell=True)
        subprocess.call("bash bash/start_blazegraph", shell=True)