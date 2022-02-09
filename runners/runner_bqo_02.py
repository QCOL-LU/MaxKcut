from max_k_cut import *
from networkx import *
import sys
import os
from read_wang_hijazi_instances import * 

num_partitions 	= int(sys.argv[1])

name 			= "Chicago" #Chicago

path 			= os.getcwd()+ "/instances/districting_2020/cities/"
final_path 		= path + name + ".txt"


graph 			= read_wang_hijazi_instance(final_path, True)

# graph 			= nx.Graph()

# edges 			= [(1,2),(1,10),(1,3),(2,10),(2,4),(2,8),(3,4),(3,9),(3,10),(4,10),(4,5),(4,6),(4,11),(5,8),(5,7),(5,11),(6,7),(6,9),(6,11),(7,11)]
# graph.add_edges_from(edges)

problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 1

problem.Params.Peel 					= False
problem.Params.Fold 					= False
problem.Params.Decompose 				= False

problem.Params.Method 					= "BQO" #A-MILO
# problem.Params.Curvature_Type 			= "convex" 			#"concave"convex
problem.Params.Curvature_Method 		= "bound"
problem.Params.Rounding_Heuristic 		= False
problem.Params.Relaxed 					= False

# problem.Params.Gurobi_Cuts 				= 3
# problem.Params.Gurobi_Presolve 			= 2
# problem.Params.Gurobi_TimeLimit			= 10
problem.Params.Gurobi_LogToConsole 		= 1
problem.Params.Gurobi_MIPGap 			= 0


problem.solve()