import networkx as nx 
from gurobipy import *
import itertools
from copy import deepcopy
from time import time
import numpy as np



#================================================================================================
# Fix twin vertices Rehfeldt et al. (2022) Faster exact solution of sparse MaxCut and QUBO problems
#================================================================================================
def twin_fix(self):

	self.preprocessing_start_time				= time()

	#--------------------------------------------------------------------------------------------
	# Find vertex pairs with distance less and equal than 2
	#--------------------------------------------------------------------------------------------
	vertex_pair_ditance 	= {(vertex1, vertex2):3 for (vertex1, vertex2) in itertools.combinations(self.vertices, 2)} 
	
	for vertex1 in self.vertices:

		ditance 			= nx.single_source_dijkstra_path_length(self.graph, vertex1, cutoff=2, weight=1)
		for (vertex2, value) in ditance.items():
			vertex_pair_ditance[(vertex1, vertex2)] 	= value
			vertex_pair_ditance[(vertex2, vertex1)] 	= value

	#--------------------------------------------------------------------------------------------
	# Initialize the variables (graph.nodes[v1][fixed-in] = v2  =>  v1 is removed and v2 is updated) fixed-in = parent
	#--------------------------------------------------------------------------------------------
	is_fixed 				= False
	fixed_graph 			= deepcopy(self.graph)

	nx.set_node_attributes(self.graph, -1, "fixed-in")


	for (vertex1, vertex2) in itertools.combinations(self.vertices, 2):

		if vertex_pair_ditance[(vertex1, vertex2)] > 2: 
			continue


		if (not fixed_graph.has_node(vertex1)) or (not fixed_graph.has_node(vertex2)): continue


		neighbors_of_vertex1 		= set(fixed_graph.neighbors(vertex1) )
		neighbors_of_vertex2 		= set(fixed_graph.neighbors(vertex2) )

		if vertex2 in neighbors_of_vertex1:
			vertices_are_linked 	= True
			neighbors_of_vertex1.remove(vertex2)
			neighbors_of_vertex2.remove(vertex1)

		else: 
			vertices_are_linked 	= False

		if not (neighbors_of_vertex1 == neighbors_of_vertex2): continue

		stop 			= False

		neighbor 		= next(iter(neighbors_of_vertex1)) 

		alpha 			= fixed_graph.edges[vertex1, neighbor]["weight"] / fixed_graph.edges[vertex2, neighbor]["weight"]

		for neighbor in neighbors_of_vertex1:
			ratio 		= fixed_graph.edges[vertex1, neighbor]["weight"] / fixed_graph.edges[vertex2, neighbor]["weight"]

			if abs(alpha - ratio) > 1e-8: 
				stop 	= True
				continue

		if stop == True: continue

		if vertices_are_linked == True:
			stop 		= (np.sign(fixed_graph.edges[vertex1, vertex2]["weight"]) * np.sign(alpha) ) > 0

		if stop == True: continue

		if np.sign(alpha) < 0 and self.num_partitions >= 3: continue

	
		#-----------------------------------------------------------------------------------------
		# Fold vertices (graph.nodes[v1][fixed-in] = v2  =>  v1 is removed and v2 is updated)
		#-----------------------------------------------------------------------------------------
		is_fixed 						= True
		self.applied_operation  		= "twin-fix"

		self.graph.nodes[vertex1]["fixed-in"] 					= vertex2

		for neighbor in neighbors_of_vertex1:
			neighbor_weight1 									= fixed_graph.edges[vertex1, neighbor]["weight"] * np.sign(alpha)
			neighbor_weight2 									= fixed_graph.edges[vertex2, neighbor]["weight"]
			neighbor_weight_new 								= neighbor_weight1 + neighbor_weight2
			fixed_graph.edges[vertex2, neighbor]["weight"] 		= neighbor_weight_new

			if abs(neighbor_weight_new) < 1e-8:
				fixed_graph.remove_node(vertex2, neighbor)

			fixed_graph.nodes[vertex2]["pos-weight"]			+= max(neighbor_weight_new, 0) - max(neighbor_weight2, 0)
			fixed_graph.nodes[vertex2]["neg-weight"]			+= min(neighbor_weight_new, 0) - min(neighbor_weight2, 0)


		if self.fixed_vertex == vertex1:
			self.fixed_vertex 									= vertex2
			partition 											= self.graph.nodes[vertex1]["partition"] 
			vertex2_partition 									= partition if np.sign(alpha) > 0 else next(iter(set(self.partitions) - set([partition]))) 
			self.graph.nodes[vertex2]["partition"] 				= vertex2_partition 
			fixed_graph.nodes[vertex2]["partition"] 			= vertex2_partition

		fixed_graph.remove_node(vertex1)

	self.preprocessing_end_time				= time()
	self.preprocessing_total_time 			= self.preprocessing_end_time  - self.preprocessing_start_time


	return (is_fixed, fixed_graph)




#================================================================================================
# Update the parent of the current node in the decomposition tree with the folding operation applied
#================================================================================================
def update_parent_twin_fix(self):

	self.parent.preprocessing_start_time			= time()
	upper_bound_change								= 0

	for vertex in self.parent.vertices:
		fixed_vertex								= vertex
		fixed_vertex_temp 							= self.parent.graph.nodes[vertex]["fixed-in"]
		partition_iter 								= 0

		while(fixed_vertex_temp != -1):
			neighbor 								= next(nx.common_neighbors(self.parent.graph, fixed_vertex_temp, fixed_vertex))
			alpha_sign 								= np.sign(self.parent.graph.edges[fixed_vertex_temp, neighbor]["weight"]) \
														* np.sign(self.parent.graph.edges[fixed_vertex, neighbor]["weight"])

			partition_iter 							+= (1 if alpha_sign < 0 else 0)

			fixed_vertex 							= fixed_vertex_temp
			fixed_vertex_temp 						= self.parent.graph.nodes[fixed_vertex_temp]["fixed-in"]

		partition 									= self.graph.nodes[fixed_vertex]["partition"]

		if (vertex != fixed_vertex) and partition_iter % 2 == 1:
			upper_bound_change 	 					+= (self.parent.graph.nodes[vertex]["pos-weight"] + self.parent.graph.nodes[vertex]["neg-weight"] )
			partition 								= next(iter(set(self.partitions) - set([partition])))
		
		self.parent.graph.nodes[vertex]["partition"]= partition

	self.parent.upper_bound 						= self.upper_bound + upper_bound_change
	
	self.parent.preprocessing_end_time				= time()
	self.parent.preprocessing_total_time 			+= self.parent.preprocessing_end_time  - self.parent.preprocessing_start_time

	


