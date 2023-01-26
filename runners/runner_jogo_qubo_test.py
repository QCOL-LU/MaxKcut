import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')

from max_k_cut import *
from networkx import *
import sys
import os
from read_hojny_et_al_instances import * 

methods = {0: "MISDO", 1: "BQO", 2: "A-MILO", 3: "RP-MILO", 4: "P-MILO", 5: "C-QUBO", 6: "CR-QUBO"}

 

name			= sys.argv[1]

ext 			= ".col"

name_ext 		= name + ext

# path 			= os.getcwd()+ "/../instances/rehfeldt/misc/" + name_ext
path 			= os.getcwd()+ "/../instances/jogo_qubo_naive/n250_d005/" + name_ext

num_partitions 	= int(sys.argv[2])
method 			= methods[int(sys.argv[3])]
penalty_setting = 0


graph 			= read_hojny_et_el_instance(path)




problem			= Instance(graph, name_specifier=name )


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 1

# problem.Params.Peel 					= False
problem.Params.Fold 					= False
# problem.Params.Decompose 				= False

# problem.Params.Twin_Fix 				= False
# problem.Params.Edge_Based_Fix 			= False
problem.Params.Rehfeldt_Fix 			= False
problem.Params.Lange_Fold 				= False

problem.Params.Method 					= method 

# problem.Params.Gurobi_LogToConsole 		= 1 


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
# problem.Params.Gurobi_TimeLimit			= 3600

problem.Params.Gurobi_TimeLimit			= 1 #600

problem.solve()

