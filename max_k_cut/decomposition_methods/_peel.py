import networkx as nx 
from copy import deepcopy
from time import time



#================================================================================================
# Peeling of the graph with positive weight iteratively
#================================================================================================
def peel(self):

	self.preprocessing_start_time				= time()
	#--------------------------------------------------------------------------------------------
	# Extend the input graph by connecting vertices with negative edge weights to auxiliary vertices
	#--------------------------------------------------------------------------------------------
	extended_graph 			= deepcopy(self.graph)


	for vertex in self.vertices:
		
		if self.graph.nodes[vertex]["neg-weight"] < 0:

			for vertex_ind in range(self.num_partitions ):
				neighbor	= "aux" + str(vertex_ind)
				extended_graph.add_edge(vertex, neighbor, weight=1)

	
			
	#--------------------------------------------------------------------------------------------
	# Find the k-core of the input graph based on extended graph
	#--------------------------------------------------------------------------------------------

	extended_graph_k_core	= set([])
	self.removed_vertices 	= []

	while True:
		temp_removed_vertices 	= [vertex for vertex, degree in extended_graph.degree() if degree < self.num_partitions]

		if not temp_removed_vertices: break

		for vertex in temp_removed_vertices:
			extended_graph.remove_node(vertex)
			if isinstance(vertex, int) or "aux" not in vertex:
				self.removed_vertices.insert(0, vertex)
			


	
	extended_graph_k_core 	= set(extended_graph.nodes)
	new_vertices 			= set(["aux" + str(vertex_ind) for vertex_ind in range(self.num_partitions )])
	self.k_core  			= list((extended_graph_k_core).difference(new_vertices))

	#--------------------------------------------------------------------------------------------
	# Determining the vertices in the k-core
	#--------------------------------------------------------------------------------------------
	nx.set_node_attributes(self.graph, False, "is_k_core")

	for vertex in self.k_core:
		self.graph.nodes[vertex]["is_k_core"] 	= True


	is_reduced 				= not (len(self.k_core) == self.num_vertices)

	self.applied_operation 	= "peel" if is_reduced == True else self.applied_operation

	
	self.preprocessing_end_time					= time()
	self.preprocessing_total_time 				= self.preprocessing_end_time  - self.preprocessing_start_time

	return is_reduced





#================================================================================================
# Update the parent of the current node in the decomposition tree with the k-core reduction operation applied
#================================================================================================
def update_parent_peel(self):

	self.parent.preprocessing_start_time	= time()
	parent_fixed_vertex_partition 			= self.parent.graph.nodes[self.parent.fixed_vertex]["partition"]

	#--------------------------------------------------------------------------------------------
	# Specify the partition of visited vertices
	#--------------------------------------------------------------------------------------------
	for vertex in self.vertices:
		self.parent.graph.nodes[vertex]["partition"]		= self.graph.nodes[vertex]["partition"]

	unvisited_vertices						= list(set(self.parent.vertices) - set(self.parent.k_core))
	visited_vertices 						= deepcopy(self.parent.k_core) 

	self.parent.upper_bound					= self.upper_bound + sum(self.parent.graph.nodes[vertex]["pos-weight"] for vertex in unvisited_vertices) - sum(self.parent.graph.edges[edge]["weight"] for edge in self.parent.graph.subgraph(unvisited_vertices).edges)


	#--------------------------------------------------------------------------------------------
	# Specify the partition of unvisited vertices
	#--------------------------------------------------------------------------------------------
	for vertex in self.parent.removed_vertices:
		neighbor_vertices 					= [neighbor for neighbor in self.parent.graph.neighbors(vertex) if neighbor in visited_vertices]
		visited_partitions 					= [self.parent.graph.nodes[vertex]["partition"] for vertex in neighbor_vertices]

		vertex_eligible_partitions 			= list(set(self.partitions).difference(set(visited_partitions) ) )


		self.parent.graph.nodes[vertex]["partition"] 	= vertex_eligible_partitions[0]

		unvisited_vertices.remove(vertex)
		visited_vertices.append(vertex)




	#--------------------------------------------------------------------------------------------
	# Change the partition of the fixed-vertex such that it becomes what it was intended to
	#--------------------------------------------------------------------------------------------		
	selected_partition 			= self.parent.graph.nodes[self.parent.fixed_vertex]["partition"]
	
	for vertex in self.parent.vertices:
		vertex_partition 		= self.parent.graph.nodes[vertex]["partition"]
		
		if vertex_partition == selected_partition:	
			self.parent.graph.nodes[vertex]["partition"] 	= parent_fixed_vertex_partition

		elif vertex_partition == parent_fixed_vertex_partition:
			self.parent.graph.nodes[vertex]["partition"] 	= selected_partition


	self.parent.preprocessing_end_time		= time()
	self.parent.preprocessing_total_time 	+= self.parent.preprocessing_end_time  - self.parent.preprocessing_start_time


	

