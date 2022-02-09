import os.path
from os import walk
import os 
import networkx as nx



def read_hojny_et_el_instance(path):

	graph 			= nx.Graph()

	with open(path) as f:
		lines 		= f.read().splitlines()

		for line in lines:
			line_list	= line.split()
			if len(line_list) > 0 and line_list[0].lower() == "e":
				weight 						= 1 if len(line_list) == 3 else line_list[-1]
				vertex1 					= int(line_list[1])
				vertex2 					= int(line_list[2])
				weight 						= float(weight)
				graph.add_edge(vertex1, vertex2, weight=weight)

	return graph


