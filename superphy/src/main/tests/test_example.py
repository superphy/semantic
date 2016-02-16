import unittest
from flask import url_for

class FlaskClientTestCase(unittest.TestCase):
    def setUp(self):
        self.assertTrue("FOO" in "FOOBAR")

    def tearDown(self):
    	self.assertTrue("BING" in "FOOBAR")

    def test_home_page(self):
        self.assertTrue("HELLO" in "HELLO WORLD")