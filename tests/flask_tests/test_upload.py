import os
import unittest

from flask import json, jsonify, url_for

from apptester import AppTester
from SuperPhy.blueprints.upload import views as module


class DataTestCase(AppTester):
    """
    tests the upload namespace
    """
    def test_something(self):
        """
        See if an upload endpoint can be called (Temporary test)
        """
        resp = self.client.get(
            url_for('upload.genome_example'),
            headers=self.get_api_headers("", ""),
            data=None
        )
        self.assertEqual(200, resp.status_code)

if __name__ == '__main__':
    unittest.main()