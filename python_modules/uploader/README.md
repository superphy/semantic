#Uploader Branch
To use with your copy of Blazegraph (in particular, the NanoSPARQL server version):

- run Blazegraph
- make a new database named superphy, with inferences (and optionally, full text index)
- run upload_all_ontologies (in ontology_uploader)

Functions are still scattered throughout the project, including queries and modifies. Use them as templates.


Proposed workflow:
- the parser takes in either raw data files (.gb, .xml) or JSON files from the existing project's parsers
    and stores information in appropriate variables
- depending on what data is available, different functions from the object_to_rdf_converter will be called
    to be stored in a triplestore graph
- the graph will be serialized and sent to the database
- a testing framework will confirm proper integration with the database, and alerts users to errors in uploading
    if it failed, and the cause

- ideally, this will be somehow integrated with a web framework


To use:

Download and run the NanoSPARQLServer from https://wiki.blazegraph.com/wiki/index.php/NanoSparqlServer
Go to Blazegraph (http://localhost:9999/bigdata/) and go to Namespaces.
Make a namespace called superphy with Modes for inference and full text index checked off.

In command line, run `pip install -r /path/to/requirements.txt` (the requirements.txt file is part of the github files

Now, go to main.py and hit run.

If you want to play with the parser and uploader, the relevant files are BioParser, object_to_rdf_converter, and
ontology_uploader
