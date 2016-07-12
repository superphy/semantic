#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest

from SuperPhy.models.upload.metadata_validation import MetadataValidator

class MetadataValidatorTestCase(unittest.TestCase):
    def setUp(self):
        self.case = MetadataValidator()

    def tearDown(self):
        pass


    def test_something(self):
        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
