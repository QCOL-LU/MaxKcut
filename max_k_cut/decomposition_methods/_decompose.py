import networkx as nx 
from copy import deepcopy
from time import time

		
#================================================================================================
# Decompose the graph based on the biconnected components
#================================================================================================
def decompose(self):

	self.preprocessing_start_time				= time()

	for vertex in self.graph.nodes:
		self.graph.nodes[vertex]["component"] 	= []
	
	self.component_iter 						= 0

	not_visited_biconnected_components 			= {ind: component for (ind, component) in enumerate(nx.biconnected_components(self.graph))}
	cut_vertices 								= set(nx.articulation_points(self.graph))
	not_visited_cut_vertices 					= deepcopy(cut_vertices)

	is_decomposed 								= False
	self.fixed_vertices 						= [self.fixed_vertex]
	self.ordered_components						= []


	if not_visited_biconnected_components:

		#----------------------------------------------------------------------------------------
		# Choose the first component based on the fixed vertex
		#----------------------------------------------------------------------------------------
		for (ind, component) in not_visited_biconnected_components.items():
			if self.fixed_vertex in component:
				visted_vertices 			= deepcopy(component)
				self.ordered_components.append(deepcopy(component))
				del not_visited_biconnected_components[ind]
				break

		
		for vertex in visted_vertices:
			self.graph.nodes[vertex]["component"].append(self.component_iter)

		
		#----------------------------------------------------------------------------------------
		# Choose the components iteratively by traversing the block tree
		#----------------------------------------------------------------------------------------
		while (len(visted_vertices) != self.num_vertices):
			intersection		= visted_vertices.intersection(not_visited_cut_vertices)
			
			if not intersection:
				key 			= next(iter(not_visited_biconnected_components))
				fixed_vertex 	= next(iter(not_visited_biconnected_components[key]))

				self.graph.nodes[fixed_vertex]["partition"] 	= 0
			else:
				fixed_vertex 	= next(iter(intersection))
				not_visited_cut_vertices.remove(fixed_vertex)


			removed_ind_components 			= []

			for (ind, component) in not_visited_biconnected_components.items():
				if fixed_vertex in component:
					
					self.component_iter 	+= 1

					self.ordered_components.append(deepcopy(component))
					self.fixed_vertices.append(fixed_vertex)
					
					for vertex in component:
						self.graph.nodes[vertex]["component"].append(self.component_iter)
						visted_vertices.add(vertex)

					removed_ind_components.append(ind)
			
			for ind in removed_ind_components:
				del not_visited_biconnected_components[ind]


	if self.component_iter > 0:
		self.applied_operation 				= "decompose"
		self.num_biconnected_component		= len(self.ordered_components)

		is_decomposed 						= True

	self.preprocessing_end_time				= time()
	self.preprocessing_total_time 			+= self.preprocessing_end_time  - self.preprocessing_start_time


	return is_decomposed



#================================================================================================
# Update the parent of the current node in the decomposition tree with the decompose operation applied
#================================================================================================
def update_parent_decompose(self):
	self.parent.preprocessing_start_time			= time()

	for vertex in self.vertices:
		self.parent.graph.nodes[vertex]["partition"]		= self.graph.nodes[vertex]["partition"]

	self.parent.preprocessing_end_time				= time()
	self.parent.preprocessing_total_time 			+= self.parent.preprocessing_end_time  - self.parent.preprocessing_start_time

	self.parent.upper_bound							+= self.upper_bound
