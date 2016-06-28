import os
import unittest

from flask import json, jsonify, url_for

from apptester import AppTester
from SuperPhy.blueprints.data import views as module



class DataTestCase(AppTester):
    """
    tests the data namespace

    You need to have your blazegraph running, and populated with test-data in
    order for these functions to work.
    """
    def test_meta(self):
        """
        General query that returns all genomes and their metadata.
        """
        resp = self.client.get(
            url_for('data.meta'),
            headers=self.get_api_headers("", ""),
            data=None
        )
        self.assertEqual(200, resp.status_code)

    def test_genomes(self):
        """
        General query that returns all genomes and their metadata.
        """
        resp = self.client.get(
            url_for('data.genomes'),
            headers=self.get_api_headers("", ""),
            data=None
        )
        self.assertEqual(200, resp.status_code)

    def test_genome(self):
        """
        @data.route('/genome/<genomeid>', methods=['GET', 'POST'])

        Returns the metadata of a particular genome in json format.
        """
        resp = self.client.get(
            url_for('data.genomes', genomeid='AAJT00000000'),
            headers=self.get_api_headers("", ""),
            data=None
        )
        self.assertEqual(200, resp.status_code)


    def test_genes(self):
        """
        @data.route('/genes', methods=['GET', 'POST'])

        General query that returns all genes and their metadata.
        """
        resp = self.client.get(
            url_for('data.genes'),
            headers=self.get_api_headers("", ""),
            data=None
        )
        self.assertEqual(200, resp.status_code)

    def test_vfs(self):
        """
        @data.route('/vf', methods=['GET', 'POST'])

        General query that returns all virulence factors.
        """
        resp = self.client.get(
            url_for('data.vfs'),
            headers=self.get_headers(),
            data=None
        )
        self.assertEqual(200, resp.status_code)

    def test_amrs(self):
        """
        @data.route('/amr', methods=['GET', 'POST'])

        General query that returns all antimicrobial resistance genes.
        """
        resp = self.client.get(
            url_for('data.amrs'),
            headers=self.get_api_headers("", ""),
            data=None
        )
        self.assertEqual(200, resp.status_code)

    def test_gene(self):
        """
        @data.route('/gene/<geneid>', methods=['GET', 'POST'])

        Returns the metadata of a particular gene in json format.
        """
        resp = self.client.get(
            url_for('data.amrs', geneid='aatA'),
            headers=self.get_api_headers("", ""),
            data=None
        )
        self.assertEqual(200, resp.status_code)

if __name__ == '__main__':
    unittest.main()