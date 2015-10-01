__author__ = 'ubiquitin'

from rdflib import Graph, Namespace, Literal

g = Graph()

n = Namespace("https://github.com/superphy/")

g.add ( (n.Ecoli, n.has_genome, Literal("O157_H7_Sakai")) )

g.serialize(destination="results.txt",format="turtle")