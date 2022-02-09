from max_k_cut import *
from networkx import *
import sys
import os
from read_hojny_et_al_instances import * 

path 			= os.getcwd()+ "/" + sys.argv[1] #instances_conncut_random
name_extension 	= sys.argv[2]

final_path 		= path 

name 			= name_extension[:-4] + "_peel_decompose_weighted"

num_partitions 	= int(name.split("_")[-1]) if "band" in name else int(sys.argv[3])


# num_partitions 	= int(sys.argv[1])

# name 			= "i160-014"
# extension 		= ".stp"

# path 			= os.getcwd()+ "/instances/hojny_et_al_instances/I160/" #instances_conncut_random
# final_path 		= path + name + extension


graph 			= read_hojny_et_el_instance(final_path)

# graph 			= nx.Graph()

# edges 			= [(1,2),(1,10),(1,3),(2,10),(2,4),(2,8),(3,4),(3,9),(3,10),(4,10),(4,5),(4,6),(4,11),(5,8),(5,7),(5,11),(6,7),(6,9),(6,11),(7,11)]
# graph.add_edges_from(edges)

problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions
# problem.Params.Verbosity 				= 1

# problem.Params.Peel 					= False
problem.Params.Fold 					= False
# problem.Params.Decompose 				= False

problem.Params.Method 					= "BQO" #A-MILO
problem.Params.Curvature_Type 			= "convex" 			#"concave"convex
problem.Params.Curvature_Method 		= "bound"
problem.Params.Rounding_Heuristic 		= False
problem.Params.Relaxed 					= False

# problem.Params.Gurobi_Cuts 				= 3
# problem.Params.Gurobi_Presolve 			= 2
problem.Params.Gurobi_TimeLimit			= 5
# problem.Params.Gurobi_LogToConsole 		= 1



problem.solve()