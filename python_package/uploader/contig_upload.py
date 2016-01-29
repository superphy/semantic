#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This module looks for eligible genomes and checks NCBI Genbank and WGS databases for relevant sequence metadata,
such as FASTA files and descriptions. Genomes are eligible if they do not already have a valid core or WGS sequence
associated with them, or if they have already been processed and found to have an invalid genome.

Classes:
	ContigUploader: Class for retrieving and uploading contig metadata of eligible genomes
	ContigsWrapper: A wrapper class that stores information for uploading and verification of contig metadata
"""

from ftplib import FTP
import gc
import gzip
import hashlib
import sys
import traceback
from urllib2 import HTTPError

from Bio import Entrez, SeqIO
from rdflib import Graph

from _sparql import check_NamedIndividual, find_missing_sequences
from _utils import generate_output, generate_path, strip_non_alphabetic
from classes import Contig
from blazegraph_upload import BlazegraphUploader
from sequence_validation import SequenceValidator

reload(sys)
sys.setdefaultencoding("utf-8")


class ContigUploader(object):
	""" Class for retrieving and uploading contig metadata of eligible genomes
	"""
	def upload_missing_contigs(self):
		""" Compiles a list of genomes with missing and unvalidated sequences and uploads them to Blazegraph.
		"""

		# find missing_sequences() returns a list of two-tuples (genome, accession number)
		for (genome, accession) in find_missing_sequences():
			contigswrapper = ContigsWrapper(genome, accession)
			try:
				self.get_seqdata(contigswrapper)
				if contigswrapper.dict["is_from"] == "PLASMID":
					self.upload(contigswrapper, self.plasmid_rdf)
				else:
					SequenceValidator(contigswrapper).validate()
					self.upload(contigswrapper, self.plasmid_rdf)
				gc.collect
			except TypeError:
				self.error_logging(contigswrapper)

	def load_contigs(self, handle, contigswrapper):
		"""
		Attempts to load contigs to Blazegraph.

		Args:
			contigswrapper: a ContigsWrapper instance for storing contig metadata
		"""
		sequencedata = SeqIO.parse(handle, 'fasta')
		contigs = [] # list of 2-tuples (accession name, seq)
		print "in load contigs outside for loop"

		for record in sequencedata:
			print "in load contigs"
			accession_name = record.name.split("|")[3].split(".")[0]
			if "complete" in record.description.lower():
				accession_name = accession_name + "_closed"

			if check_NamedIndividual(accession_name):
				print "%s already in Blazegraph." % accession_name
				raise TypeError
			else:
				print "Genome complete in load contigs"
				contigs.append((accession_name, str(record.seq)))

				if "plasmid" in record.description.lower():
					contigswrapper.dict["is_from"] = 'PLASMID'

				contigswrapper.add_contigs(contigs)


	def get_seqdata(self, contigswrapper):
		"""
		Args:
			contigswrapper: a ContigsWrapper instance that holds contig metadata for a genome

		Returns: a BLAST record for self.load_contigs to use
		"""
		Entrez.email = "superphy.info@gmail.com"
		handle = None
		i = 0

		while(i<3):
			try:
				handle = Entrez.efetch(db="nuccore", id=contigswrapper.genome, rettype="fasta", retmode="text")
				for record in SeqIO.parse(handle, 'fasta'):
					if "complete" in record.description.lower():
						contigswrapper.dict["is_from"] = "CORE"
						handle = Entrez.efetch(db="nuccore", id=contigswrapper.genome, rettype="fasta", retmode="text")
						self.load_contigs(handle, contigswrapper)
						break
					else:
						self.download_file(strip_non_alphabetic(str(contigswrapper.genome)), 'fsa_nt.gz')
						with open(generate_path('tmp/loading.fasta'), 'rb') as handle:
							contigswrapper.dict["is_from"] = "WGS"
							self.load_contigs(handle, contigswrapper)
			except HTTPError:
				i += 1
				print i
				continue
			break
		try:
			handle is None
		except NameError:
			raise TypeError("Could not retrieve file for analysis")


	def upload(self, contigswrapper, func):
		"""Uploads contig data to Blazegraph based on inputted function argument

		Args: 
			contigswrapper: ContigsWrapper instance that is a wrapper for contig metadata
			func: function that handles RDF conversion based on the type of genome
		"""
		g = Graph()
		for (accession_name, seq) in contigswrapper.contigs:
			contig_rdf = Contig(g, **contigswrapper.generate_kwargs(accession_name, seq))
			func(contigswrapper, contig_rdf)

			BlazegraphUploader().upload_data(generate_output(g))


	## FIX THESE AFTER LUNCH
	def plasmid_rdf(self, contigswrapper, contig_rdf):
		"""Sets up RDF triples for a plasmid sequence and its metadata. Marks genome as possessing a valid sequence if
		the accession id matches the genome id (i.e. metadata validation has not merged genomes with the same biosample
		id together yet)

		Args:
			contigswrapper: a ContigsWrapper instance storing contig data that would otherwise be a data clump
			contig_rdf: an initalized RDF converter for Contig data (i.e. a Contig instance)
		"""
		contig_rdf.rdf()

		if contigswrapper.accession == contigswrapper.genome:
			contig_rdf.add_seq_validation(True)

	def nonplasmid_rdf(self, contigswrapper, contig_rdf):
		"""Sets up RDF triples for a nonplasmid sequence and its metadata in accordance with its sequence validation
		 results

		Args:
			contigswrapper: a ContigsWrapper instance storing contig data that would otherwise be a data clump
			contig_rdf: an initalized RDF converter for Contig data (i.e. a Contig instance)
		"""
		contig_rdf.add_seq_validation(contigswrapper.valid)

		if contigswrapper.valid:
			contig_rdf.rdf()
			contig_rdf.add_hits(contigswrapper.hits)


	def error_logging(self, contigswrapper):
		"""Logs errors regarding contig uploading to a file, for manual curation.

		Args:
			contigswrapper: a ContigsWrapper instance storing sequence-related data that would otherwise be a data clump
		"""
		with open(generate_path("outputs/seq_errors.txt"), "a") as f:
			f.write("Genome: %s - Accession: %s.\n" % (contigswrapper.genome, contigswrapper.accession))
			f.write("%s \n ================================ \n\n" % traceback.format_exc())
		print "%s - %s: The records for this sequence are not retrievable." % (contigswrapper.genome, contigswrapper.accession)


	def download_file(self, id, filetype):
		"""Downloads the gzip file with the correct id and filetype and unzips it and transfers its contents into a
		temporary FASTA file for further processing. If no files on the server match, returns a TypeError.

		Args:
			id(str): a WGS project ID, composed of only alphabetics
			filetype(str): the type of file to be found. 'fsa_nt.gz' is the default, but there are other options for
				amino acids and other formats
		"""
		ftp = FTP('bio-mirror.jp.apan.net')
		ftp.login('anonymous', 'superphy.info@gmail.com')
		ftp.cwd('pub/biomirror/genbank/wgs')

		filenames = ftp.nlst()
		filename = [s for s in filenames if id in s and filetype in s]

		if len(filename) is not 1:
			raise TypeError("No files could be found for download.")
		else:
			ftp.retrbinary('RETR ' + filename[0], open(generate_path('tmp/loading.gz'), 'wb').write)
			with gzip.open(generate_path('tmp/loading.gz')) as fasta, \
					open(generate_path('tmp/loading.fasta'), 'wb') as output:
				output.write(fasta.read())


class ContigsWrapper(object):
	"""A wrapper class that holds metadata for contigs of a certain genome, in addition 
	"""
	def __init__(self, genome, accession):
		"""Initializes the class with the necessary fields. The dict is used to construct the kwargs to pass into
		classes.Contigs for uploading into Blazegraph

		Args:
			genome(str): the accession of the genome that the sequence is associated with
			accession(str): the accession that the sequence is identified by
		"""
		self.name = str(accession) + "_seq"
		self.accession = str(accession)
		self.genome = str(genome)
		self.hits = None
		self.valid = None
		self.bp = None
		self.contigs = None
		self.numcontigs = None
		self.checksum = None


		keys = ["name", "genome", "sequence", "is_from"]
		self.dict = {key: None for key in keys}
		self.dict["genome"] = self.genome

	def add_contigs(self, contigs):
		"""
		Args:
			contigs: a list of 2-tuples where 2-tuples are in the form (accession_name(str), seq(str))
		"""
		self.contigs = contigs
		self.numcontigs = len(contigs)
		self.bp = sum(len(seq) for (accession_name, seq) in contigs)
		self.generate_checksum()

	def generate_checksum(self):
		"""Generates an MD5 checksum from the sorted contigs in the sequence (sorted from smallest to largest; order
		matters).
		"""

		sorted_contigs = sorted(self.contigs, key=lambda contig: len(contig[1]))
		seqhash = hashlib.md5()
		for (accession_name, seq) in sorted_contigs:
			seqhash.update(str(seq))
		self.checksum = seqhash.hexdigest()

	def generate_kwargs(self, accession_name, seq):
		"""Returns a dict of kwargs to pass into the classes.Contig constructor. If any dict entries are NoneType,
		raise an error indicating the first found that is missing.
		"""
		self.dict["name"] = accession_name
		self.dict["sequence"] = seq
		for key, value in self.dict.iteritems():
			if not value:
				raise TypeError("Missing contig metadata: %s" % key)

		return self.dict


if __name__ == "__main__":
	#import pdb; pdb.set_trace()
	ContigUploader().upload_missing_contigs()


