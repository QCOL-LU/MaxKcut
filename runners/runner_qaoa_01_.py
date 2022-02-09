from max_k_cut import *
from networkx import *
import sys

seed 			= int(sys.argv[1])
num_vertices 	= int(sys.argv[2])

name 			= "seed" + str(seed)  
graph 			= nx.Graph()

# # edges 			= [(1,2),(1,10),(1,3),(2,10),(2,4),(2,8),(3,4),(3,9),(3,10),(4,10),(4,5),(4,6),(4,11),(5,8),(5,7),(5,11),(6,7),(6,9),(6,11),(7,11)]
# # graph.add_edges_from(edges)


# edges 			= [(1,2,1), (1,3,1), (2,3,2), (2,5,1), (2,4,1),(5,4,-2.5), (3,5,1)]
edges 			= [(1,2,-1), (1,3,-1), (2,3,3), (2,4,1), (3,4,1)]
# edges 			= [(1,2,-1), (1,3,-1), (2,3,-1), (2,4,3), (3,4,-1), (1,4,-1)]
# edges 			= [(1,2,-1), (1,3,-1), (1,4,-1), (1,5,-1), (2,3,-1), (2,4,-1), (2,5,-1), (3,4,-1), (3,5,-1), (4,5,-1)]

graph.add_weighted_edges_from(edges)
graph 			= nx.gnp_random_graph(num_vertices, 0.8, seed=seed)
problem			= Instance(graph, name_specifier=name)


problem.Params.Num_Partitions 			= 2
problem.Params.K_Core 					= False
problem.Params.Clique_Constraints		= False
problem.Params.Verbosity 				= 2

problem.Params.Fold 					= False
problem.Params.Decompose 				= False
problem.Params.QAOA_Verbosity 			= 2
# problem.Params.QAOA_Scipy_Optimizer 	= "shgo"

problem.Params.Method 					= "BQO"
# problem.Params.Curvature_Type 			= "convex"
# problem.Params.Curvature_Type 		= "concave"
problem.Params.Curvature_Method 		= "power_iteration"
problem.Params.Rounding_Heuristic 		= False
# problem.Params.Symmetry_Breaking 		= True

problem.solve()