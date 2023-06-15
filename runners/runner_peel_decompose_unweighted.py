from max_k_cut import *
from networkx import *
import sys
import os
from read_wang_hijazi_instances import * 


path 			= os.getcwd()+ "/"+sys.argv[1] #instances_conncut_random
name_extension 	= sys.argv[2]

final_path 		= path 
num_partitions 	= int(sys.argv[3])

name 			= name_extension[:-4] + "_peel_decompose_unweighted"





graph 			= read_wang_hijazi_instance(final_path, False)

problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions
# problem.Params.Peel 					= False
problem.Params.Verbosity 				= 1

problem.Params.Fold 					= False
# problem.Params.Decompose 				= False

problem.Params.Method 					= "BQO" #A-MILO
problem.Params.Curvature_Type 			= "convex" 			#"concave"convex
problem.Params.Curvature_Method 		= "numpy"

problem.Params.Rounding_Heuristic 		= False
problem.Params.Relaxed 					= False

# problem.Params.Gurobi_Cuts 				= 3
problem.Params.Gurobi_Presolve 			= 2
problem.Params.Gurobi_TimeLimit			= 10



problem.solve()

