from max_k_cut import *
from networkx import *
import sys
import os
from read_wang_hijazi_instances import * 


path 			= os.getcwd()+ "/"+sys.argv[1] #instances_conncut_random
name_extension 	= sys.argv[2]

final_path 		= path 

name 			= name_extension[:-4] 

num_partitions 	= int(name.split("_")[-1]) if "band" in name else int(sys.argv[3])

name 			= name + "_con"

name 			= name_extension[:-4] + "_weighted_con"

# name 			= "random10_100" #Chicago

# path 			= os.getcwd()+ "/instances/wang_hijazi_instances/not_in_paper/"
# # path 			= os.getcwd()+ "/instances/districting_2020/cities/"

# final_path 		= path + name + ".txt"

# name 			= name + "_weighted_con_preprocess_th10" #_weighted_con_dec_peel
# num_partitions 	= int(sys.argv[1])


graph 			= read_wang_hijazi_instance(final_path, False)

problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 1

# problem.Params.Peel 					= False
# problem.Params.Fold 					= False
# problem.Params.Decompose 				= False

problem.Params.Method 					= "BQO" #A-MILO
problem.Params.Curvature_Type 			= "convex" 			#"concave"convex
problem.Params.Curvature_Method 		= "bound"

problem.Params.Rounding_Heuristic 		= False
problem.Params.Relaxed 					= False

# problem.Params.Gurobi_Cuts 				= 3
# problem.Params.Gurobi_Presolve 			= 2
# problem.Params.Gurobi_TimeLimit			= 10 #3600
# problem.Params.Gurobi_Threads			= 10
# problem.Params.Gurobi_LogToConsole 		= 1


problem.solve()

