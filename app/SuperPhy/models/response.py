"""
models.py
"""
#pylint: disable=C0103, W0406

from collections import defaultdict
from flask import Flask, abort, request, jsonify, g, url_for, current_app
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.httpauth import HTTPBasicAuth
from passlib.apps import custom_app_context as pwd_context
from itsdangerous import (TimedJSONWebSignatureSerializer
                          as Serializer, BadSignature, SignatureExpired)

from flask import send_file
import StringIO

class Response(object):
    """
    This is a helper for formatting the sparql responses for the user.

    """
    @classmethod
    def csv(cls, dictionary):
        """
        Converts a Blazegraph response dictionary into a tab separated values.

        Because we want our first key to be Accession id number, make sure your
        sparql query is formed such that that is the first key.
        """
        response = dictionary

        #the dict stores all the metadata in the heads,vars list
        metadata = response["head"]["vars"]

        #we want to print every genome as a new row, and the metadata in the same order
        #using the accession number as the genome name
        csv = "\t".join(metadata) + "\n"
        #the JSON stores every genome as an entry in results,bindings
        genomeObj = response["results"]["bindings"]

        for genome in genomeObj:
            #start the next line of output
            line = []

            #get all of the other metadata for printing
            for m in metadata:
                #some values have include a newline
                #we would prefer not to have that for our table output
                nextM = genome[m]["value"].replace("\n", " ")
                line.append(nextM)

            #end the current line, prepare for the next one
            print line
            csv += "\t".join(line) + "\n"

        return csv

    @classmethod
    def csvfile(cls, dictionary):
        """
        Returns response as a file.
        """
        data = Response.csv(dictionary)
        strIO = StringIO.StringIO()
        strIO.write(str(data))
        strIO.seek(0)
        return send_file(strIO,
                         attachment_filename="metadata.txt",
                         as_attachment=True)

    @classmethod
    def default(cls, results, extra=None):
        """
        This is a template. Copy this over to another function if you are
        making a different response object.
        """
        response = {}

        response = results

        if extra is not None:
            response.update(extra)
        return jsonify(response)

    @classmethod
    def bulk_download(cls, results, extra=None):
        """
        """
        response = {}
        response['headers'] = results['head']['vars']
        response['rows'] = []
        bindings = results['results']['bindings']
        for binding in bindings:
            row = {}
            for item in response['headers']:
                if binding.has_key(item):
                    row[item] = binding[item]['value']
                else:
                    row[item] = ""
            response['rows'].append(row)

        if extra is not None:
            response.update(extra)
        return jsonify(response)

    @classmethod
    def shorten(cls, results, extra=None):
        """
        This is the default response from sparql. This turns the layered
        response into a shorter dictionary object
        Example

        input:
            {'head': {'vars': ['s', 'p', 'o']}, 'results': {'bindings': [{'p': {'value': 'FOOP'}, 's': {'value': 'FOOS'}, 'o': {'value': 'FOOO'}}, {'p': {'value': 'FOOP2'}, 's': {'value': 'FOOS2'}, 'o': {'value': 'FOOO2'}}]}}
        output:
        {'vars': ['s', 'p', 'o'], 'results': [{'p': 'FOOP', 's': 'FOOS', 'o': 'FOOO'}, {'p': 'FOOP2', 's': 'FOOS2', 'o': 'FOOO2'}]}
        """
        response = {}

        response['vars'] = results.get('head').get('vars')
        bindings = results['results']['bindings']
        response['rows'] = []
        for key in bindings:
            row = {}
            for item in response['vars']:
                row[item] = key[item]['value']
            response['rows'].append(row)

        if extra is not None:
            response.update(extra)
        return jsonify(response)

    @classmethod
    def format_gene_search(cls, results, genomeDict, extra=None):
        """
        Formats the response from the gene search query for the front-end in the form of a dictionary:
        {genome1: {gene1: count, gene2: count}, genome2: {gene1: count...} ...}

        Extra arguments:
            genomeDict: a preformed dictionary based on selected genes and genomes.
        """
        response = {}
        response = results

        bindings = results['results']['bindings']
        for binding in bindings:
            accession = binding['Genome']['value'].split("#")[1]
            gene_name = binding['Gene_Name']['value']
            try:
                genomeDict[accession][gene_name] += 1
            except KeyError:
                "Genome or gene doesn't exist in dictionary"

        if extra is not None:
            response.update(extra)

        return jsonify(genomeDict)

    @classmethod
    def format_categories(cls, results, extra=None):
        """
        Formats the response for categories to use for the front-end.
        """
        response = {}
        response = results

        categoryDict = defaultdict(list)

        bindings = results['results']['bindings']
        for binding in bindings:
            category = binding["Category"]['value']
            subcategory = binding["Subcategory"]['value']
            if subcategory not in categoryDict[category]:
                categoryDict[category].append(subcategory) 
        return jsonify(categoryDict)
