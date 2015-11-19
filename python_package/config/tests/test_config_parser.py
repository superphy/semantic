#!/usr/bin/env python

"""Tests for config.parser module

Example:
        $ python test_config_parser.py

"""

import unittest
from superphy.config import parser

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "CC-SY"
__version__ = "4.0"
__maintainer__ = "Matt"
__email__ = "matthew.whiteside@canada.ca"


class ConfigParserTestCase(unittest.TestCase):

    # def shell_source(self, script):
    #     """
    #     Emulate the action of "source" in bash,
    #     settings some environment variables.

    #     """
    #     pipe = subprocess.Popen(["/bin/bash","-c",". %s; env" % script], stdout=subprocess.PIPE, shell=True)
    #     output = pipe.communicate()[0]
    #     print output
    #     env = dict((line.split("=", 1) for line in output.splitlines()))

    #     environ.update(env)


    def setUp(self):
        pass
       
    def tearDown(self):
        pass
       

    def test_parse_superphy_environment_variables(self):
        """
        Enviroment variables are initizialized with values
        matching INI-config file settings. Environment initizialization
        happens in the postactivate.sh script called in venv/bin/activate.

        This test run parser.read() to check that the parser works and that
        the environment contains the proper variables.

        """
        try:
            props = parser.read()
            self.assertIsNotNone(props['rdf_url'], "Environment variable SUPERPHY_RDF_URL found")

        except parser.SuperphyConfigError, e:
            self.fail("parser.read() raised SuperphyConfigError! {}".format(str(e)))



if __name__ == '__main__':
    unittest.main()
