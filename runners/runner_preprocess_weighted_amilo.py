import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')

from max_k_cut import *
from networkx import *
import sys
import os
from read_wang_hijazi_instances import * 


# path 			= os.getcwd()+ "/" + sys.argv[1] #instances_conncut_random

path 			= "/home/ramin/MaxKcut/instances/wang_hijazi_instances/band/" + "band250_4.txt"
name_extension 	= ".txt"#sys.argv[2]

final_path 		= path 

name 			= name_extension[:-4] 

num_partitions 	= int(name.split("_")[-1]) if "band" in name else 4 # int(sys.argv[3])

# name 			= name + "_weighted"

# name 			= "Chicago" #Chicago

# path 			= os.getcwd()+ "/instances/districting_2020/cities/"
# final_path 		= path + name + ".txt"

# name 			= name + "_weighted"
# num_partitions 	= 2 #int(sys.argv[1])


graph 			= read_wang_hijazi_instance(final_path, True)

problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= num_partitions

problem.Params.Verbosity 				= 1

# problem.Params.Peel 					= False
# problem.Params.Fold 					= False
# problem.Params.Decompose 				= False

problem.Params.Method 					= "A-MILO" #A-MILO BQO
# problem.Params.Curvature_Type 			= "convex" 			#"concave"convex
# problem.Params.Curvature_Method 		= "bound"
problem.Params.Gurobi_LogToConsole 		= 1 

problem.Params.Rounding_Heuristic 		= False
problem.Params.Relaxed 					= False
problem.Params.Clique_Constraints 		= True

# problem.Params.Gurobi_Cuts 				= 3
# problem.Params.Gurobi_Presolve 			= 2
problem.Params.Gurobi_TimeLimit			= 3600



problem.solve()

