import networkx as nx 
from gurobipy import *
import itertools
from copy import deepcopy
from time import time
import numpy as np



#================================================================================================
# Fixed vertices based on Lange et al. (2019) Combinatorial Persistency Criteria for Multicut and Max-Cut.
#================================================================================================
def edge_based_fix(self):

	self.preprocessing_start_time				= time()

	#--------------------------------------------------------------------------------------------
	# Initialize the variables (graph.nodes[v1][fixed-in] = v2  =>  v1 is removed and v2 is updated) fixed-in = parent
	#--------------------------------------------------------------------------------------------
	is_fixed 				= False
	fixed_graph 			= deepcopy(self.graph)

	nx.set_node_attributes(self.graph, -1, "fixed-in")


	for (vertex1, vertex2) in self.edges:


		if (not fixed_graph.has_edge(vertex1, vertex2)) : continue


		edge_weight 				= fixed_graph.edges[(vertex1, vertex2)]["weight"]
		total_weight_vertex1 		= fixed_graph.nodes[vertex1]["pos-weight"] - fixed_graph.nodes[vertex1]["neg-weight"] - abs(edge_weight)
		total_weight_vertex2 		= fixed_graph.nodes[vertex2]["pos-weight"] - fixed_graph.nodes[vertex2]["neg-weight"] - abs(edge_weight)

		edge_fixed 					= abs(edge_weight) >= min(total_weight_vertex1, total_weight_vertex2)

		if not edge_fixed: continue

		if (edge_weight > 0 and self.num_partitions > 2): continue

	
		#-----------------------------------------------------------------------------------------
		# Fold vertices (graph.nodes[v1][fixed-in] = v2  =>  v1 is removed and v2 is updated)
		#-----------------------------------------------------------------------------------------
		is_fixed 					= True
		self.applied_operation  	= "edge-based-fix"

		self.graph.nodes[vertex1]["fixed-in"] 					= vertex2

		neighbors_of_vertex1 		= set(fixed_graph.neighbors(vertex1) )
		neighbors_of_vertex1.remove(vertex2)

		for neighbor in neighbors_of_vertex1:
			neighbor_weight1 									= fixed_graph.edges[vertex1, neighbor]["weight"] * (-np.sign(edge_weight))
			if not fixed_graph.has_edge(vertex2, neighbor):
				fixed_graph.add_edge(vertex2, neighbor, weight=0)

			neighbor_weight2 									= fixed_graph.edges[vertex2, neighbor]["weight"] 
			neighbor_weight_new 								= neighbor_weight1 + neighbor_weight2
			fixed_graph.edges[vertex2, neighbor]["weight"] 		= neighbor_weight_new

			if abs(neighbor_weight_new) < 1e-8:
				fixed_graph.remove_edge(vertex2, neighbor)

			fixed_graph.nodes[vertex2]["pos-weight"]			+= max(neighbor_weight_new, 0) - max(neighbor_weight2, 0)
			fixed_graph.nodes[vertex2]["neg-weight"]			+= min(neighbor_weight_new, 0) - min(neighbor_weight2, 0)


		if self.fixed_vertex == vertex1:
			self.fixed_vertex 									= vertex2
			partition 											= self.graph.nodes[vertex1]["partition"] 
			vertex2_partition 									= partition if edge_weight < 0 else next(iter(set(self.partitions) - set([partition]))) 
			self.graph.nodes[vertex2]["partition"] 				= vertex2_partition 
			fixed_graph.nodes[vertex2]["partition"] 			= vertex2_partition

		fixed_graph.remove_node(vertex1)

	self.preprocessing_end_time				= time()
	self.preprocessing_total_time 			= self.preprocessing_end_time  - self.preprocessing_start_time


	return (is_fixed, fixed_graph)




#================================================================================================
# Update the parent of the current node in the decomposition tree with the folding operation applied
#================================================================================================
def update_parent_edge_based_fix(self):

	self.parent.preprocessing_start_time			= time()

	upper_bound_change								= 0

	for vertex in self.parent.vertices:
		fixed_vertex								= vertex
		fixed_vertex_temp 							= self.parent.graph.nodes[vertex]["fixed-in"]
		partition_iter 								= 0


		while(fixed_vertex_temp != -1):
			edge_weight 							= self.parent.graph.edges[fixed_vertex_temp, fixed_vertex]["weight"]
			partition_iter 							+= (1 if edge_weight > 0 else 0)

			
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

	


