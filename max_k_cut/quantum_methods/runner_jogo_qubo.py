import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')

from max_k_cut import *
from networkx import *
import sys
import os
from read_hojny_et_al_instances import * 

methods = {0: "MISDO", 1: "BQO", 2: "A-MILO", 3: "RP-MILO", 4: "P-MILO", 5: "C-QUBO", 6: "CR-QUBO"}

path 			= sys.argv[1] 

name_ext		= sys.argv[2]
name, ext 		= name_ext.rsplit('.', 1)

num_partitions 	= int(sys.argv[3])
method 			= methods[int(sys.argv[4])]
penalty_setting = int(sys.argv[5])




final_path 		= path #+ name_ext

graph 			= read_hojny_et_el_instance(final_path)


if penalty_setting == 0:
	name 		+= "_tight" 
elif penalty_setting == 1:
	name 		+= "_adjusted_tight"
elif penalty_setting == 2:
	name 		+= "_naive"

problem			= Instance(graph, name_specifier=name )


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 1

problem.Params.Peel 					= False
problem.Params.Fold 					= False
problem.Params.Decompose 				= False
problem.Params.Fold_Twin 				= False

problem.Params.Method 					= method 

problem.Params.Gurobi_LogToConsole 		= 1 


if penalty_setting == 0:
	problem.Params.Naive				= False
	problem.Params.Adjusted_Penalty		= False

elif penalty_setting == 1:
	problem.Params.Naive				= False
	problem.Params.Adjusted_Penalty		= True

elif penalty_setting == 2:
	problem.Params.Naive				= True
	problem.Params.Adjusted_Penalty		= False


problem.Params.Rounding_Heuristic 		= False
problem.Params.Relaxed 					= False
problem.Params.Clique_Constraints 		= False

# problem.Params.Gurobi_Cuts 				= 3
# problem.Params.Gurobi_Presolve 			= 2
problem.Params.Gurobi_TimeLimit			= 3600



problem.solve()

