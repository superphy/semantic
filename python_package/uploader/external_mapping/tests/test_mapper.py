#!/usr/bin/env python

"""Tests for uploader.external_mapper.mapper module

Example:
        $ python test_mapper.py

"""

import unittest
import logging
import shutil, tempfile
from os import path
import json
from superphy.uploader.external_mapping.mapper import Mapper, SuperphyMapperError

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "ASF"
__version__ = "2.0"
__maintainer__ = "Matt"
__email__ = "matthew.whiteside@canada.ca"


class MapperTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Create tmp directory & logging 

        """
        cls.test_dir = tempfile.mkdtemp()
        print "Tmpdir: %s"%cls.test_dir

        logger = logging.getLogger('mapper.py test')
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(cls.test_dir + '/mapper_tests.log')
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        cls.logger = logger


    @classmethod
    def tearDownClass(self):
        """
        Remove tmp directory

        """
    
        #shutil.rmtree(self.test_dir)


    def setUp(self):
        self.logger.info("\n\n\n**********START OF NEW TEST*********\n\n")


    # Test initialization of Mapper object with various decision_tree JSON objects
    def test_invalid_json_in_decision_tree(self):

        # Decision tree json file
        filenm1 = path.join(self.test_dir, 'invalid_decision_tree.json')
        f1 = open(filenm1, 'w')

        # Invalid JSON string
        f1.write("{'Nope: 'hello'")
        f1.flush()
        f1.close()

        # Overrides json file
        overrides = {}
        filenm2 = path.join(self.test_dir, 'overrides.json')
        with open(filenm2, 'w') as outfile:
            json.dump(overrides, outfile)

        # Check for exception
        with self.assertRaises(SuperphyMapperError):
            Mapper(filenm1, filenm2, self.logger)


    def test_unrecognized_validation_method_in_decision_tree(self):

        # Decision tree json file
        filenm1 = path.join(self.test_dir, 'invalid_decision_tree.json')

        with open(filenm1, 'w') as f1:
            dt = { "specific_host": {
                "keep": True,
                "validation_routines": ["cant_possibly_be_a_validation_method_noway_anytime_everrrrr"],
                "cleanup_routines": ["fix_human"]
            }}
            json.dump(dt, f1)

        # Overrides json file
        overrides = {}
        filenm2 = path.join(self.test_dir, 'overrides.json')
        with open(filenm2, 'w') as outfile:
            json.dump(overrides, outfile)

        # Check for exception
        with self.assertRaises(SuperphyMapperError):
            Mapper(filenm1, filenm2, self.logger)


    def test_unrecognized_cleanup_method_in_decision_tree(self):

        # Decision tree json file
        filenm1 = path.join(self.test_dir, 'invalid_decision_tree.json')

        with open(filenm1, 'w') as f1:
            dt = { "specific_host": {
                "keep": True,
                "validation_routines": ["hosts"],
                "cleanup_routines": ["cant_possibly_be_a_cleanup_method_noway_anytime_everrrrr"]
            }}
            json.dump(dt, f1)

        # Overrides json file
        overrides = {}
        filenm2 = path.join(self.test_dir, 'overrides.json')
        with open(filenm2, 'w') as outfile:
            json.dump(overrides, outfile)

        # Check for exception
        with self.assertRaises(SuperphyMapperError):
            Mapper(filenm1, filenm2, self.logger)

    # Test duplicate hosts
    def test_multiple_hosts(self):

        # Decision tree json file
        filenm1 = path.join(self.test_dir, 'invalid_decision_tree.json')
        with open(filenm1, 'w') as f1:
            dt = { "test_host1": {
                    "keep": True,
                    "validation_routines": ["hosts"],
                    "cleanup_routines": ["fix_human"]
                    },
                "test_host2": {
                    "keep": True,
                    "validation_routines": ["hosts"],
                    "cleanup_routines": ["fix_cow"]
                    },
                }
            json.dump(dt, f1)

        # Overrides json file
        overrides = {}
        filenm2 = path.join(self.test_dir, 'overrides.json')
        with open(filenm2, 'w') as f2:
            json.dump(overrides, f2)

        # Input json file
        filenm3 = path.join(self.test_dir, 'inputs.json')
        with open(filenm3, 'w') as f3:
            meta = { "accession1": {
                    "test_host1": "human",
                    "test_host2": "cow"
                }
            }
            json.dump(meta, f3)

        # Initialize mapper
        mapper = Mapper(filenm1, filenm2, self.logger)

        # Check for exception during run
        filenm4 = path.join(self.test_dir, 'output.json')
        with self.assertRaises(SuperphyMapperError):
            mapper.run(filenm3, filenm4)
            

    # Test incompatible category
    def test_incompatible_source(self):

        # Decision tree json file
        filenm1 = path.join(self.test_dir, 'invalid_decision_tree.json')
        with open(filenm1, 'w') as f1:
            dt = { "test_host": {
                    "keep": True,
                    "validation_routines": ["hosts"],
                    "cleanup_routines": ["fix_human"]
                    },
                "test_source": {
                    "keep": True,
                    "validation_routines": ["sources"],
                    "cleanup_routines": []
                    },
                }
            json.dump(dt, f1)

        # Overrides json file
        overrides = {}
        filenm2 = path.join(self.test_dir, 'overrides.json')
        with open(filenm2, 'w') as f2:
            json.dump(overrides, f2)

        # Input json file
        filenm3 = path.join(self.test_dir, 'inputs.json')
        with open(filenm3, 'w') as f3:
            meta = { "accession1": {
                    "test_host": "human",
                    "test_source": "soil"
                }
            }
            json.dump(meta, f3)

        # Initialize mapper
        mapper = Mapper(filenm1, filenm2, self.logger)

        # Check for exception during run
        filenm4 = path.join(self.test_dir, 'output.json')
        with self.assertRaises(SuperphyMapperError):
            mapper.run(filenm3, filenm4)
            


    def test_mapper(self):

        # Use test files in etc/ directory
        etcdir = path.dirname(path.realpath(__file__)) + "/../etc/"

        dtfile = etcdir + "test_decision_tree.json"
        orfile = etcdir + "test_overrides.json"
        infile = etcdir + "test_input.json"
        outfile = etcdir + "test_output.json" # This file will get overwritten

        mapper = Mapper(dtfile, orfile, self.logger)
        mapper.run(infile, outfile)

        # Expected genome record
        expected = {
            "ANWA00000000": {
                "isolation_host": 'hsapiens'
            }
        }
        # Read output
        with open(outfile) as f:
            data = json.load(f)

            # Check output genome record has proper metadata assignments
            self.assertEqual(expected, data)




if __name__ == '__main__':
    unittest.main()