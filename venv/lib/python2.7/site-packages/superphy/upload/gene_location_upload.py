#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""This module uploads gene location information on genomes.

Classes:
	GeneLocationUploader: parses blast data results and stores info for uploading to Blazegraph.
"""

from collections import defaultdict
import sys
import traceback
import subprocess
import json
import re
import os

from rdflib import Graph
from Bio.Blast import NCBIXML

from _sparql import check_NamedIndividual, has_ref_gene, _sparql_query
from _utils import generate_output, generate_path
from classes import GeneLocation
from blazegraph_upload import BlazegraphUploader
from contig_upload import ContigUploader

__author__ = "Clarice Ng"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Clarice Ng"
__email__ = "c32ng@uwaterloo.ca"


class GeneLocationUploader(object):
	"""A class for parsing and uploading gene locations from a Blasted XML file.

	"""

	def __init__(self):
		self.dict = {}
		

	def upload(self, filename):
		"""
		"""
		self.ncbixml_parse(filename)
		print "\n BEGIN GET REFERENCES"
		reference_genes = self.get_reference_genes()
		print "\n BEGIN CREATE FASTA"
		self.create_fasta(reference_genes, "tmp/ref_sequences.fasta")
		self.create_db("samples/sample.fasta")
		self.blastn_commandline("genome_db")
		print "\n BEGIN PARSING LOCAL BLAST RESULTS"
		self.parse_result()


	def accession_name(self, contig, desc):
		"""
		Returns a string of the accession number with "_closed" appended at the end if it is a complete genome.

		Args:
			contig(str): an accession number for a contig 
			desc(str): a description of a particular contig/genome.
		"""
		name = contig

		if "complete genome" in desc:
			contig = contig + "_closed"

		return contig


	def add_contig(self, gene_name, contig):
		"""
		Adds contig to self.dict to count the number of times a gene has appeared on a certain contig.

		Args:
			gene_name(str): name of the gene
			contig(str): accession name of the contig
		"""
		if (gene_name in self.dict) and (contig in self.dict[gene_name]):
			self.dict[gene_name][contig] += 1
		else:
			self.dict[gene_name] = {}
			num_copies = self.get_num_gene_copies(gene_name, contig)
			self.dict[gene_name][contig] = num_copies


	def get_gene_name(self, s):
		"""
		Returns the gene name from a string.
		Args:
			str(s): a string (most likely from a blast query or fasta file)
		"""
		gene_name = ""
		try:
			if "ARO" in s:
				gene_name = s.split()[1].split(".")[0]
			elif "VFO" in s:
				gene_name = s.split("|")[0]
			else:
				raise ValueError

			if "/" in gene_name:
				return gene_name.split("/")[0]
			elif "(" in gene_name:
				return gene_name.split("(")[0]
			else:
				return gene_name
		except(ValueError):
			print "Not a valid AMR or virulence factor file"

		# try:
		# 	matchObj = re.match(r'(([a-z][a-zA-Z][a-z](_?[a-zA-Z][a-zA-Z]+[0-9]?|[a-zA-Z][0-9]*(_[0-9])?|[0-9]*[A-Z]*)(-[0-9])?)|z[0-9][0-9][0-9][0-9]|eCs[0-9]+)(.|\|)', s)
		# 	if matchObj:
		# 		print "match --> matchObj.group() : ", matchObj.group(0)
		# 		return matchObj.group(0)
		# 	else:
		# 		raise ValueError
		# except(ValueError):
		# 	searchObj = re.search(r'(([a-z][a-zA-Z][a-z](_?[a-zA-Z][a-zA-Z]+[0-9]?|[a-zA-Z][0-9]*(_[0-9])?|[0-9]*[A-Z]*)(-[0-9])?)|z[0-9][0-9][0-9][0-9]|eCs[0-9]+)(\.|\|)+', s)
		# 	if searchObj:
		# 		print "search --> searchObj.group() : ", searchObj.group(0)
		# 		return searchObj.group(0)
		# 	else:
		# 		print "Nothing found."

	def get_num_gene_copies(self, gene, contig):
		"""
		Queries the database to return the number of occurrences of a gene in a contig there are.

		Args:
			gene(str): the name of the gene
			contig(str): the name of the contig
		"""

		results = _sparql_query(
			'PREFIX : <https://github.com/superphy#>\n'
			'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
			'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
			'PREFIX faldo: <http://biohackathon.org/resource/faldo#>\n'
			'SELECT ?GeneLocation \n'
        	'WHERE { :%s :has_copy ?GeneLocation . '
        	':%s :has_gene ?GeneLocation . }' % (gene, contig)
    	)
		copies = 0

		for result in results["results"]["bindings"]:
			copies += 1

		return copies 


	def ncbixml_parse(self, filename):
		"""
		Takes the XML virulence factor blastn result and adds metadata to the GeneMetadata instance. 
		Also calls method to upload sequence information for genomes matched that are not already in the Blazegraph

		"""
		metadata = None
		nonuploaded_genomes = []
		missing_alignments = []

		with open(generate_path(filename)) as result_handle:
			blast_records = NCBIXML.parse(result_handle)

			E_VALUE_THRESH = 0.04
			count = 0

			for blast_record in blast_records:
				# gene_name = blast_record.query.split("|")[0]
				# if "/" in gene_name: #if gene has multiple names
				# 	gene_name = gene_name.split("/")[0]
				gene_name = self.get_gene_name(blast_record.query)

				max_percentage = -1
				min_gaps = 1000

				print gene_name
				print "count:", count

				location_name, contig_name, begin, end, sequence, ref_gene, uploaded = [None, None, None, None, None, None, None]

				for alignment in blast_record.alignments:
					alignment_descr = alignment.title.split("|")[6]
					contig_accession = alignment.title.split("|")[5].split(".")[0]
					contig_name = self.accession_name(contig_accession, alignment_descr)

					for hsp in alignment.hsps:
						percentage = hsp.score/hsp.identities * 100

						if hsp.expect < E_VALUE_THRESH and hsp.gaps <= min_gaps and percentage >= max_percentage:
							max_percentage = percentage
							min_gaps = hsp.gaps

							location_name = gene_name + "_" + contig_name + "_" + "0"

							if (not check_NamedIndividual(contig_name)) and ((contig_accession, contig_accession) not in nonuploaded_genomes):
								nonuploaded_genomes.append((contig_accession, contig_accession))
							
							begin = hsp.sbjct_start
							end = str(int(hsp.sbjct_start) + int(hsp.score) - 1)
							ref_gene = False

							if not has_ref_gene(gene_name):
								ref_gene = True

							if self.is_complete_genome(alignment_descr):
								self.add_contig(gene_name, contig_name)
								location_name = gene_name + "_" + contig_name + "_" + str(self.dict[gene_name][contig_name])

								print "Uploading", location_name
								self.create_gene_location(location_name, gene_name, contig_name, begin, end, hsp.sbjct, ref_gene)
								uploaded = True
								break

					if uploaded:
						break

				if not uploaded:
					if not location_name:
						missing_alignments.append(gene_name)
					else:
						self.add_contig(gene_name, contig_name)
						self.create_gene_location(location_name, gene_name, contig_name, begin, end, hsp.sbjct, ref_gene)

				count += 1

		with open("missing_alignments.txt", "w") as f:
			for n in missing_alignments:
				f.write("%s\n" % n)

		#ContigUploader().upload_missing_contigs(nonuploaded_genomes)

	def is_complete_genome(self, descr):
		"""
		Determines whether the given description describes a complete genomes.
		Returns a boolean.

		Args:
			descr(str): a description of a certain contig or genome
		"""
		if "complete genome" in descr:
			return True
		else:
			return False


	def add_to_graph(self, metadata):
		"""Attempts to upload data to Blazegraph via conversion to RDF and turtle. Tracks errors made during this process
		as it likely has to do with either missing curated data, a blocked or inaccessible NCBI site, and issues with
		formatting and information availability.

		Args:
			metadata (GeneMetadata): a GeneMetadata object used to store relevant key-value pairs from the parse
		"""
		if metadata:
			try:
				self.progress += 1
				print "%d: downloading files" % self.progress
				self.create_gene(metadata)
			except Exception as e:
				self.error_logging(metadata.name)


	def create_gene_location(self, name, gene, contig, begin, end, seq, ref_gene):
		""" Creates a GeneLocation object to export the data in turtle format with the appropriate RDF tags and
		uploads it to Blazegraph.

		Args:
			name (str): gene name  + "_" + contig name + "_" + occurence number
            gene (str): the name of the gene
            contig (str): name of the contig its found in
            begin (str): the beginning position of the gene in contig
            end (str): the end position of the gene in contig
            seq (str): the sequence of the gene on the contig
            is_ref_gene (boolean): signifies whether this is the reference gene location for a particular gene
		"""

		if self.check_gene_copy(gene, contig, begin, end):
			print "This copy of %s in %s is already in Blazegraph." % (gene, contig)
		else:
			g = Graph()
			GeneLocation(g, name, gene, contig, begin, end, seq, ref_gene).rdf()
			BlazegraphUploader().upload_data(generate_output(g))


	def check_gene_copy(self, gene, contig, begin, end):
		"""
		Checks to see if a certain gene copy on a contig has already been uploaded.
		"""

		results = _sparql_query(
			'PREFIX : <https://github.com/superphy#>\n'
			'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
			'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
			'PREFIX faldo: <http://biohackathon.org/resource/faldo#>\n'
			'ASK { ?begin faldo:position "%s"^^xsd:string . ?end faldo:position "%s"^^xsd:string . ?genecopy faldo:begin ?begin . ?genecopy faldo:end ?end . ?genecopy :is_gene_of :%s . :%s :has_copy ?genecopy . }\n' % (begin, end, contig, gene)
		)

		return results["boolean"]


	def get_reference_genes(self):
		"""
		Returns a list of all the reference gene instances and their sequences for analysis (from Blazegraph)
		"""

		results = _sparql_query(
			'PREFIX : <https://github.com/superphy#>\n'
			'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
			'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
			'PREFIX faldo: <http://biohackathon.org/resource/faldo#>\n'
			'SELECT ?s ?gene ?seq WHERE { ?s rdf:type faldo:Region . ?s rdf:type :reference_gene . ?s :has_sequence ?seq . ?gene :has_copy ?s .}'
		)

		return ((result["s"]["value"].rsplit("#", 1)[1], result["gene"]["value"].split("#")[1], result["seq"]["value"])
				for result in results["results"]["bindings"])


	def create_db(self, filename):
		"""
		Creates a nucleotide database from self.filename for comparing against the reference genes in Blazegraph.
		"""

		db_path = "../../../../../../blast/ncbi-blast*/bin/makeblastdb"

		subprocess.call("%s -dbtype nucl -title genome_db -out genome_db -in %s" 
						% (db_path, generate_path(filename)), shell=True)


	def create_fasta(self, genes, out_file):
		"""Writes a FASTA sequence to a file for use by the command line version of BLAST. Obtains nucleotide data from
		the sequence data object used to initialize the validator and writes each entry as a separate FASTA object.
		"""
		with open(generate_path(out_file), "w") as f:
			for (name, gene_name, seq) in genes:
				f.write(">%s\n%s\n" %(gene_name, seq))


	def blastn_commandline(self, db):
		"""
		Runs the command line blastn on the reference gene sequences using the supplied database.
		Outputs the results to tmp/results.xml

		Args:
			db(str): A file path to a database used for blast comparison.
		"""
		blastn_path = generate_path("../../../../../../blast/ncbi-blast*/bin/blastn")

		command = blastn_path
		fasta = generate_path("tmp/ref_sequences.fasta")
		path_to_db = generate_path(db)
		print command
		results = generate_path("tmp/results.xml")

		subprocess.call("%s -query %s -db %s -outfmt 5 -out %s -best_hit_score_edge 0.05 -best_hit_overhang 0.1"
						 % (command, fasta, path_to_db, results),
						shell=True)


	def parse_result(self):
		"""
		Takes the XML result from the virulence factor blast and adds metadata to the GeneMetadata instance. 
		Also calls method to upload sequence information for genomes matched that are not already in Blazegraph

		"""

		with open(generate_path("tmp/results.xml")) as result_handle:
			blast_records = NCBIXML.parse(result_handle)

			E_VALUE_THRESH = 0.04

			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < E_VALUE_THRESH:
							gene_name = blast_record.query

							contig_accession = alignment.title.split("|")[5].split(".")[0]
							contig_name = self.accession_name(contig_accession, alignment.title.split("|")[6])

							self.add_contig(gene_name, contig_name)

							name = gene_name + "_" + contig_name + "_" + str(self.dict[gene_name][contig_name])
							begin = hsp.sbjct_start
							end = str(int(hsp.sbjct_start) + int(hsp.score) - 1)
							ref_gene = False

							self.create_gene_location(name, gene_name, contig_name, begin, end, hsp.sbjct, ref_gene)


###### For Testing purposes ######

if __name__ == "__main__":
 	# For gene testing
	gmd1 = GeneLocationUploader()
	gmd1.upload('data/superphy_amr.xml')
	#gmd1.upload('data/superphy_vf.xml')