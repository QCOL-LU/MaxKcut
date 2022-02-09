import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')

from max_k_cut import *
from networkx import *
import sys
import numpy as np
import random

seed 				= int(sys.argv[1])
num_vertices 		= int(sys.argv[2])
num_partitions 		= int(sys.argv[3])
graph_density		= float(sys.argv[4])
neg_edge_percentage	= float(sys.argv[5])
penalty_increase 	= float(sys.argv[6])

noise 				= 0.01

is_connected 		= False
while (not is_connected):
	
	random.seed(seed) 
	graph 			= nx.erdos_renyi_graph(num_vertices, graph_density, seed=seed, directed=False)
	seed 			= seed + 1 
	is_connected 	= nx.is_connected(graph)


weighted_graph	= nx.Graph()
edges 			= [(edge[0], edge[1], 1 if random.random() > neg_edge_percentage else -1) for edge in graph.edges()] 


name 				= "seed" + str(seed - 1) + "_rand_p" + str(graph_density) + "_neg" + str(neg_edge_percentage) + ("naive" if penalty_increase > 1 else "tight") + "_noise_" + str(noise)
# name 				= "seed" + str(seed - 1) + "_rand_p" + str(graph_density) + "_neg" + str(neg_edge_percentage) 

weighted_graph.add_weighted_edges_from(edges)
problem			= Instance(weighted_graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 2

problem.Params.Peel 					= False
problem.Params.Fold 					= False
problem.Params.Decompose 				= False

problem.Params.Method 					= "QUBO" #R-QUBO BQO
problem.Params.Naive 					= (True if penalty_increase > 1 else False)


problem.Params.Gurobi_TimeLimit			= 3600
# problem.Params.QAOA_Optimize 			= False

problem.Params.QAOA_Scipy_Optimizer 	= "brute" #Nelder-Mead brute
# problem.Params.Penalty_Increase 		= penalty_increase

# problem.Params.QAOA_Angles				= [0.26, 0.38] #[0.26, 0.38]
problem.Params.QAOA_Num_Levels			= 1
problem.Params.QAOA_Num_Shots			= 10000
problem.Params.QAOA_Brute_Num_Samples 	= 25
problem.Params.QAOA_Gates_Noise			= noise

problem.solve()