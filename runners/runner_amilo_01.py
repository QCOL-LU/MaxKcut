from max_k_cut import *
from networkx import *
import sys

seed 			= int(sys.argv[1])
num_vertices 	= int(sys.argv[2])
num_partitions 	= int(sys.argv[3])

name 			= "seed" + str(seed)  
graph 			= nx.Graph()

graph 			= nx.gnp_random_graph(num_vertices, 0.5, seed=seed)
problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions
problem.Params.K_Core 					= False

problem.Params.Clique_Constraints		= True
problem.Params.Bicycle_Wheel_Constraints= False
problem.Params.Wheel_Constraints 		= True

problem.Params.Verbosity 				= 2

problem.Params.Fold 					= False
problem.Params.Decompose 				= False


problem.Params.Method 					= "A-MILO"
problem.Params.Gurobi_TimeLimit			= 5 * 60 

problem.solve()


