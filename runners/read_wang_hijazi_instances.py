import os.path
from os import walk
import os 
import networkx as nx
import numpy as np



def read_wang_hijazi_instance(path, is_weighted):

	graph 			= nx.Graph()
	np.random.seed(0)

	with open(path) as f:
		lines 		= f.read().splitlines()[1:]
		for line in lines:
			(vertex1, vertex2, weight) 	= line.split( )
			vertex1 					= int(vertex1)
			vertex2 					= int(vertex2)
			rand_num 					= np.random.random_sample()# - 0.5
			weight 						= float(weight) if is_weighted else np.sign(rand_num)
			graph.add_edge(vertex1, vertex2, weight=weight)


	return graph


