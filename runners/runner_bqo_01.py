from max_k_cut import *
from networkx import nx
import sys
import numpy as np
from read_wang_hijazi_instances import *

# seed 			= int(sys.argv[1])
# num_vertices 	= int(sys.argv[2])
num_partitions 	= int(sys.argv[1])
# penalty_increase= float(sys.argv[4])

# name 			= "band150_3" #Chicago

# path 			= os.getcwd()+ "/instances/wang_hijazi_instances/band/"
# final_path 		= path + name + ".txt"


graph 			= read_wang_hijazi_instance(final_path, True)

# name 			= "seed" + str(seed) + "_pen"+  str(penalty_increase) 
# graph 			= nx.Graph()

# edges 			= [(1,2),(1,10),(1,3),(2,10),(2,4),(2,8),(3,4),(3,9),(3,10),(4,10),(4,5),(4,6),(4,11),(5,8),(5,7),(5,11),(6,7),(6,9),(6,11),(7,11)]
# graph.add_edges_from(edges)


# edges 			= [(1,2,-1), (1,3,-2), (2,3,4), (2,4,1), (2,5,1), (3,4,1), (3,5,1), (4,5,1)]	#QUBO example
# edges 			= [(1,2,1), (1,3,1), (1,4,0), (2,3,-1), (2,5,1), (3,4,1), (2,4,1), (3,5,1)]		#R-QUBO example
# graph.add_weighted_edges_from(edges)


# np.random.seed(0) 
# edges = [(i, j, 2*np.random.randint(0, 2) - 1) for i in np.arange(1, num_vertices, 1) for j in np.arange(i + 1, min(i+num_partitions+1, num_vertices) + 1, 1)]

# # print(edges)
# graph.add_weighted_edges_from(edges)

print("is_chordal:", nx.is_chordal(graph))
print(hi)

# graph 			= nx.gnp_random_graph(num_vertices, 0.8, seed=seed)
problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 2

problem.Params.Peel 					= False
# problem.Params.Fold 					= False
problem.Params.Decompose 				= False

problem.Params.Method 					= "R-QUBO" #R-QUBO BQO

# problem.Params.Curvature_Type 			= "convex" 			#"concave"
# problem.Params.Curvature_Method 		= "power_iteration"

problem.Params.Gurobi_TimeLimit			= 3600
# problem.Params.Gurobi_LogToConsole 		= 1

problem.Params.QAOA_Scipy_Optimizer 	= "brute" #Nelder-Mead brute
problem.Params.Penalty_Increase 		= penalty_increase

# problem.Params.QAOA_Angles				= [2.79, 3.49, 1.40, 5.59]#[1.74141958, 1.69713506, 0.56961772, 1.5603471,  0.63246351, 0.97906295] #[np.pi/4 for i in range(6)]
# problem.Params.QAOA_Angles				= [0.50, 0.60, 2.07, 4.80]
problem.Params.QAOA_Num_Levels			= 1
# problem.Params.QAOA_Num_Shots			= 500


problem.solve()