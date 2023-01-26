import networkx as nx 
from gurobipy import *
import itertools
from copy import deepcopy
from time import time
import numpy as np



#================================================================================================
# Fold the graph based on Lange et al. (2019) Combinatorial Persistency Criteria for Multicut and Max-Cut.
#================================================================================================
def lange_fold(self):

	self.preprocessing_start_time				= time()

	#--------------------------------------------------------------------------------------------
	# Initialize the variables (graph.nodes[v1][folded-in] = v2  =>  v1 is removed and v2 is updated) folded-in = parent
	#--------------------------------------------------------------------------------------------
	is_folded 				= False
	folded_graph 			= deepcopy(self.graph)
	temp_graph 				= deepcopy(self.graph)

	nx.set_node_attributes(self.graph, -1, "folded-in")


	for vertex1 in self.vertices:

		if vertex1 not in list(temp_graph.nodes): continue

		neighbors_of_vertex1 		= list(temp_graph.neighbors(vertex1) )

		neighbor_subraph 			= temp_graph.subgraph(neighbors_of_vertex1)

		is_vertex1_removed 			= False		

		for (vertex2, vertex3) in list(neighbor_subraph.edges):

			triangle 				= [vertex1, vertex2, vertex3]

			for (tri_vertex1, tri_vertex2) in itertools.combinations(triangle, 2):
				if  (not folded_graph.has_edge(tri_vertex1, tri_vertex2)): continue

				tri_vertex3 		= next(iter(set(triangle) - set([tri_vertex1, tri_vertex2])))

				if  ((not folded_graph.has_edge(tri_vertex1, tri_vertex3)) 
					or (not folded_graph.has_edge(tri_vertex2, tri_vertex3)) ): continue 

				weight_edge12 		= max(folded_graph.edges[tri_vertex1, tri_vertex2]["weight"], 0) 
				weight_edge13 		= max(folded_graph.edges[tri_vertex1, tri_vertex3]["weight"], 0)
				weight_edge23 		= max(folded_graph.edges[tri_vertex2, tri_vertex3]["weight"], 0)

				right_hand_side1 	= folded_graph.nodes[tri_vertex1]["pos-weight"] - folded_graph.nodes[tri_vertex1]["neg-weight"] 
				right_hand_side2 	= folded_graph.nodes[tri_vertex2]["pos-weight"] - folded_graph.nodes[tri_vertex2]["neg-weight"] 


				condition1 			= 2 * (weight_edge12 + weight_edge13) >= right_hand_side1
				condition2 			= 2 * (weight_edge12 + weight_edge23) >= right_hand_side2

				
				if condition1 and condition2:

					#----------------------------------------------------------------------------
					# Fold vertices (graph.nodes[v1][folded-in] = v2  =>  v1 is removed and v2 is updated)
					#----------------------------------------------------------------------------
					is_folded 						= True
					self.applied_operation  		= "lange-fold"

					self.graph.nodes[tri_vertex1]["folded-in"] 				= tri_vertex2

					is_vertex1_removed 				= (vertex1 == tri_vertex1)

					print(tri_vertex1, tri_vertex2, tri_vertex3)

					neighbors_tri_vertex1 			= set(folded_graph.neighbors(tri_vertex1)) - set([tri_vertex2])

					for neighbor in neighbors_tri_vertex1:

						neighbor_weight1 									= folded_graph.edges[tri_vertex1, neighbor]["weight"] 
						
						if not folded_graph.has_edge(tri_vertex2, neighbor):
							folded_graph.add_edge(tri_vertex2, neighbor, weight=0)

						neighbor_weight2 									= folded_graph.edges[tri_vertex2, neighbor]["weight"]
						neighbor_weight_new 								= neighbor_weight1 + neighbor_weight2
						folded_graph.edges[tri_vertex2, neighbor]["weight"] = neighbor_weight_new

						if abs(neighbor_weight_new) < 1e-8:
							folded_graph.remove_edge(vertex2, neighbor)

						folded_graph.nodes[tri_vertex2]["pos-weight"]		+= max(neighbor_weight_new, 0) - max(neighbor_weight2, 0)
						folded_graph.nodes[tri_vertex2]["neg-weight"]		+= min(neighbor_weight_new, 0) - min(neighbor_weight2, 0)


					if self.fixed_vertex == tri_vertex1:
						self.fixed_vertex 									= tri_vertex2
						partition 											= self.graph.nodes[tri_vertex1]["partition"] 
						self.graph.nodes[tri_vertex2]["partition"] 			= partition 
						folded_graph.nodes[tri_vertex2]["partition"] 		= partition

					folded_graph.remove_node(tri_vertex1)
					temp_graph.remove_node(tri_vertex1)

					break
				
			if is_vertex1_removed: break
		
		if vertex1 in set(temp_graph.nodes):
			temp_graph.remove_node(vertex1)


	self.preprocessing_end_time				= time()
	self.preprocessing_total_time 			= self.preprocessing_end_time  - self.preprocessing_start_time


	return (is_folded, folded_graph)




#================================================================================================
# Update the parent of the current node in the decomposition tree with the folding operation applied
#================================================================================================
def update_parent_lange_fold(self):

	self.parent.preprocessing_start_time			= time()

	for vertex in self.parent.vertices:
		folded_vertex								= vertex
		folded_vertex_temp 							= self.parent.graph.nodes[vertex]["folded-in"]

		while(folded_vertex_temp != -1):
			folded_vertex 							= folded_vertex_temp
			folded_vertex_temp 						= self.parent.graph.nodes[folded_vertex_temp]["folded-in"]

		partition 									= self.graph.nodes[folded_vertex]["partition"]
		
		self.parent.graph.nodes[vertex]["partition"]= partition

	self.parent.upper_bound 						= self.upper_bound
	
	self.parent.preprocessing_end_time				= time()
	self.parent.preprocessing_total_time 			+= self.parent.preprocessing_end_time  - self.parent.preprocessing_start_time

	


