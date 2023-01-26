import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')

from max_k_cut import *
from networkx import *
import sys
import os
from read_hojny_et_al_instances import * 

methods = {0: "MISDO", 1: "BQO", 2: "A-MILO", 3: "RP-MILO", 4: "P-MILO"}

instance_name 	= "miles750.col"

path 			= "../instances/conv1/n100_d015/" + instance_name # sys.argv[1] # 

name_ext		=  instance_name # sys.argv[2] #
name, ext 		= name_ext.rsplit('.', 1) 

num_partitions 	= int(sys.argv[3])
method 			= methods[int(sys.argv[4])]


final_path 		= path #+ name_ext

graph 			= read_hojny_et_el_instance(final_path)


problem			= Instance(graph, name_specifier=name + "_concave")


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 1

problem.Params.Peel 					= False
problem.Params.Fold 					= False
problem.Params.Decompose 				= False
# 
problem.Params.Method 					= method #A-MILO BQO
problem.Params.Curvature_Type 			= "concave" 			#"concave"convex
problem.Params.Curvature_Method 		= "mosek"
problem.Params.Gurobi_LogToConsole 		= 1 

problem.Params.Rounding_Heuristic 		= False
problem.Params.Relaxed 					= True
problem.Params.Clique_Constraints 		= False

# problem.Params.Gurobi_Cuts 				= 3
# problem.Params.Gurobi_Presolve 			= 2
problem.Params.Gurobi_TimeLimit			= 3600



problem.solve()

