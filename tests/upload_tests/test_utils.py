#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import inspect
import os
import unittest

from SuperPhy.models.upload import _utils

class utilsTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_generate_path(self):
        (frame, filepath, line_number, function_name, lines, index) = inspect.stack()[0]
        expected_dir = os.path.dirname(filepath)
        generated = _utils.generate_path("asdf")
        generated_dir = os.path.dirname(generated)

        self.assertEqual(expected_dir,generated_dir)

        filename = os.path.basename(generated)

        self.assertEqual(filename, "asdf")

    def test_strip_non_alphabetic(self):
        test_string = "abcdefghi"
        result = _utils.strip_non_alphabetic(test_string)
        self.assertEqual(test_string, result)

        test_string = "abcd890"
        result = _utils.strip_non_alphabetic(test_string)
        self.assertEqual("abcd",result)

        test_string = "123456789"
        result = _utils.strip_non_alphabetic(test_string)
        self.assertEqual("", result)

    def test_strip_non_numeric(self):
        test_string = "123456789"
        result = _utils.strip_non_numeric(test_string)
        self.assertEqual(test_string, result)

        test_string = "abcd890"
        result = _utils.strip_non_numeric(test_string)
        self.assertEqual("890",result)

        test_string = "abcdefghi"
        result = _utils.strip_non_numeric(test_string)
        self.assertEqual("", result)

    def test_generate_output(self):
        pass

    def test_generate_file_output(self):
        pass

if __name__ == '__main__':
    unittest.main()
