import networkx as nx 
from gurobipy import *
import itertools
from copy import deepcopy
from time import time
import numpy as np



#================================================================================================
# Fix an edge Rehfeldt et al. (2022) Faster exact solution of sparse MaxCut and QUBO problems
#================================================================================================
def rehfeldt_fix(self):

	self.preprocessing_start_time				= time()

	#--------------------------------------------------------------------------------------------
	# Initialize the variables (graph.nodes[v1][fixed-in] = v2  =>  v1 is removed and v2 is updated) fixed-in = parent
	#--------------------------------------------------------------------------------------------
	is_fixed 				= False
	fixed_graph 			= deepcopy(self.graph)
	temp_graph 				= deepcopy(self.graph)

	nx.set_node_attributes(self.graph, -1, "fixed-in")


	for vertex1 in self.vertices:

		neighbors_of_vertex1 		= list(fixed_graph.neighbors(vertex1) )

		neighbor_subraph 			= temp_graph.subgraph(neighbors_of_vertex1)		

		is_vertex1_removed 			= False

		for (vertex2, vertex3) in neighbor_subraph.edges:

			triangle 				= [vertex1, vertex2, vertex3]

			for (tri_vertex1, tri_vertex2) in itertools.permutations(triangle, 2):

				if (not fixed_graph.has_edge(tri_vertex1, tri_vertex2)): continue

				tri_vertex3 		= next(iter(set(triangle) - set([tri_vertex1, tri_vertex2])))

				if  ((not fixed_graph.has_edge(tri_vertex1, tri_vertex3)) 
					or (not fixed_graph.has_edge(tri_vertex2, tri_vertex3)) ): continue 

				weight_edge12 		= max(fixed_graph.edges[tri_vertex1, tri_vertex2]["weight"], 0) 
				weight_edge13 		= max(fixed_graph.edges[tri_vertex1, tri_vertex3]["weight"], 0)
				weight_edge23 		= min(fixed_graph.edges[tri_vertex2, tri_vertex3]["weight"], 0)


				if not (min(weight_edge12, weight_edge13) > 0 > weight_edge23): continue

				right_hand_side1 	= fixed_graph.nodes[tri_vertex1]["pos-weight"] - fixed_graph.nodes[tri_vertex1]["neg-weight"] 
				right_hand_side2 	= fixed_graph.nodes[tri_vertex2]["pos-weight"] - fixed_graph.nodes[tri_vertex2]["neg-weight"] 


				condition1 			= 2 * (weight_edge12 + weight_edge13) >= right_hand_side1
				condition2 			= 2 * (weight_edge12 - weight_edge23) >= right_hand_side2

				
				if condition1 and condition2:

					#----------------------------------------------------------------------------
					# Fold vertices (graph.nodes[v1][fixed-in] = v2  =>  v1 is removed and v2 is updated)
					#----------------------------------------------------------------------------
					is_fixed 						= True
					self.applied_operation  		= "rehfeldt-fix"
					
					is_vertex1_removed 				= (vertex1 == tri_vertex1)

					neighbors_tri_vertex1 			= set(fixed_graph.neighbors(tri_vertex1)) - set([tri_vertex2])

					self.graph.nodes[tri_vertex1]["fixed-in"] 				= tri_vertex2


					for neighbor in neighbors_tri_vertex1:
						neighbor_weight1 									= -fixed_graph.edges[tri_vertex1, neighbor]["weight"]

						if not fixed_graph.has_edge(tri_vertex2, neighbor):
							fixed_graph.add_edge(tri_vertex2, neighbor, weight=0)

						neighbor_weight2 									= fixed_graph.edges[tri_vertex2, neighbor]["weight"]
						neighbor_weight_new 								= neighbor_weight1 + neighbor_weight2
						fixed_graph.edges[tri_vertex2, neighbor]["weight"]  = neighbor_weight_new

						if abs(neighbor_weight_new) < 1e-8:
							fixed_graph.remove_edge(vertex2, neighbor)

						fixed_graph.nodes[tri_vertex2]["pos-weight"]		+= max(neighbor_weight_new, 0) - max(neighbor_weight2, 0)
						fixed_graph.nodes[tri_vertex2]["neg-weight"]		+= min(neighbor_weight_new, 0) - min(neighbor_weight2, 0)


					if self.fixed_vertex == tri_vertex1:
						self.fixed_vertex 									= tri_vertex2
						partition 											= self.graph.nodes[vertex1]["partition"] 
						tri_vertex2_partition 								= next(iter(set(self.partitions) - set([partition]))) 
						self.graph.nodes[tri_vertex2]["partition"] 			= tri_vertex2_partition 
						fixed_graph.nodes[tri_vertex2]["partition"] 		= tri_vertex2_partition

					fixed_graph.remove_node(tri_vertex1)
					temp_graph.remove_node(tri_vertex1)


					break
			if is_vertex1_removed: break
		
		if vertex1 in set(temp_graph.nodes):
			temp_graph.remove_node(vertex1)


	self.preprocessing_end_time				= time()
	self.preprocessing_total_time 			= self.preprocessing_end_time  - self.preprocessing_start_time


	return (is_fixed, fixed_graph)






#================================================================================================
# Update the parent of the current node in the decomposition tree with the folding operation applied
#================================================================================================
def update_parent_rehfeldt_fix(self):

	self.parent.preprocessing_start_time			= time()

	upper_bound_change								= 0

	for vertex in self.parent.vertices:
		fixed_vertex								= vertex
		fixed_vertex_temp 							= self.parent.graph.nodes[vertex]["fixed-in"]

		while(fixed_vertex_temp != -1):
			fixed_vertex 							= fixed_vertex_temp
			fixed_vertex_temp 						= self.parent.graph.nodes[fixed_vertex_temp]["fixed-in"]
			edge_weight 							= self.parent.graph.edges[fixed_vertex_temp, fixed_vertex]["weight"]

		partition 									= self.graph.nodes[fixed_vertex]["partition"]
		edge_weight 								= 0
		
		if fixed_vertex != vertex:
			edge_weight 	= self.parent.graph.edges[vertex, fixed_vertex]["weight"]
			neighbor 		= next(nx.common_neighbors(self.parent.graph, vertex, fixed_vertex))
			alpha_sign 		= np.sign(self.parent.graph.edges[vertex, neighbor]["weight"]) * np.sign(self.parent.graph.edges[fixed_vertex, neighbor]["weight"])
			partition 		= partition if alpha_sign > 0 else next(iter(set(self.partitions) - set([partition])))
		
		self.parent.graph.nodes[vertex]["partition"]= partition
		upper_bound_change 	+= min(edge_weight, 0)

	self.parent.upper_bound 						= self.upper_bound
	
	self.parent.preprocessing_end_time				= time()
	self.parent.preprocessing_total_time 			+= self.parent.preprocessing_end_time  - self.parent.preprocessing_start_time

	


