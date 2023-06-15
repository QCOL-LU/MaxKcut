import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')

from max_k_cut import *
from networkx import *
import sys
import os
from read_hojny_et_al_instances import * 

methods = {0: "MISDO", 1: "BQO", 2: "A-MILO", 3: "RP-MILO", 4: "P-MILO", 5: "MISDO2"}

# instance_name 	= "queen8_12.col"

path 			=  sys.argv[1] # "../instances/or_letters_instances/n050_d030/" + instance_name #

name_ext		=  sys.argv[2] # instance_name # 
name, ext 		= name_ext.rsplit('.', 1) 

num_partitions 	= int(sys.argv[3])
method 			= methods[int(sys.argv[4])]


final_path 		= path #+ name_ext

graph 			= read_hojny_et_el_instance(final_path)


# problem			= Instance(graph, name_specifier=name + "_concave")
problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 1

problem.Params.Peel 					= False
problem.Params.Fold 					= False
problem.Params.Decompose 				= False
# 
problem.Params.Method 					= method #A-MILO BQO
# problem.Params.Curvature_Type 			= "concave" 			#"concave"convex
# problem.Params.Curvature_Method 		= "mosek"
# problem.Params.Gurobi_LogToConsole 		= 1 

# problem.Params.Rounding_Heuristic 		= False
# problem.Params.Relaxed 					= True
# problem.Params.Clique_Constraints 		= False

# problem.Params.Gurobi_Cuts 				= 3
# problem.Params.Gurobi_Presolve 			= 2
problem.Params.Gurobi_TimeLimit			= 3600



problem.solve()

