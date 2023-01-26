import networkx as nx 
from gurobipy import *
import itertools
from copy import deepcopy
from time import time
import numpy as np

# import matplotlib.pyplot as plt

def folded_subgraph_solver(self, graph, folded_vertices, add_const=False):

	vertex1, vertex2 	= folded_vertices
	print(vertex1, vertex2, list(graph.nodes))

	with Env(empty=True) as env:
		env.setParam('LogToConsole', 0)
		# sys.stdout.write(2*"\033[F\033[K")
		env.start()
		with Model(env=env) as model:

			x 		= model.addVars(list(graph.nodes), self.partitions, vtype=GRB.BINARY, name="x")

			model.setObjective(	quicksum(graph.edges[edge]["weight"] * x[edge[0],partition]*x[edge[1],partition] for edge in graph.edges for partition in self.partitions) \
									, GRB.MINIMIZE)	


			#-----------------------------------------------------------------------------------
			# Constraints
			#-----------------------------------------------------------------------------------
			model.addConstrs(quicksum(x[vertex, partition] for partition in self.partitions) == 1 for vertex in graph.nodes)

			if add_const:
				
				model.addConstrs(x[vertex1, partition] == x[vertex2, partition] for partition in self.partitions)

			model.Params.NonConvex 			= 2
			model.Params.Threads			= self.Params.Gurobi_Threads

			model.optimize()

			is_the_same 					= True

			obj_value 						= model.objVal

			if not add_const:
				is_the_same 				= all([ abs(x[vertex1, partition].x  - x[vertex2, partition].x) <= 1e-4 for partition in self.partitions])

	return (obj_value, is_the_same)				


#================================================================================================
# Fold the graph
#================================================================================================
def fold(self):

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
	# Initialize the variables (graph.nodes[v1][folded-in] = v2  =>  v1 is removed and v2 is updated) folded-in = parent
	#--------------------------------------------------------------------------------------------
	is_folded 				= False
	folded_graph 			= deepcopy(self.graph)

	nx.set_node_attributes(self.graph, -1, "folded-in")


	for (vertex1, vertex2) in itertools.combinations(self.vertices, 2):

		if vertex_pair_ditance[(vertex1, vertex2)] > 2: continue


		if (not folded_graph.has_node(vertex1)) or (not folded_graph.has_node(vertex2)): continue

		common_neighbors 	 		= list( nx.common_neighbors(folded_graph, vertex1, vertex2))		

		# extended_common_neighbors 	= common_neighbors + [vertex1, vertex2]

		# subgraph_edges				= list(folded_graph.subgraph(extended_common_neighbors).edges)


		boundary_edges				= list(nx.edge_boundary(folded_graph, [vertex1, vertex2]))

		# all_edges 					= subgraph_edges + boundary_edges
		
		new_graph 					= folded_graph.edge_subgraph(boundary_edges)

		# nx.draw(new_graph, with_labels=True)
		# new_name 					= self.name + "_" + str(vertex1) + "_" + str(vertex2) 
		# plt.savefig(self.figure_path + "/"+ new_name + ".png", dpi=300)
		# plt.clf()

	

		obj_value, can_be_folded 	= self.folded_subgraph_solver(new_graph, (vertex1, vertex2))
		# print(vertex1, vertex2, obj_value, can_be_folded )


		if not can_be_folded:
			new_obj_value, is_the_same 		= self.folded_subgraph_solver(new_graph, (vertex1, vertex2), True)

			can_be_folded			= new_obj_value <= obj_value	

			# print(new_obj_value, can_be_folded )



		

		#-----------------------------------------------------------------------------------------
		# Fold vertices (graph.nodes[v1][folded-in] = v2  =>  v1 is removed and v2 is updated)
		#-----------------------------------------------------------------------------------------
		if can_be_folded == True:
			print(vertex1, vertex2)
			is_folded 						= True
			self.applied_operation  		= "fold"

			self.graph.nodes[vertex1]["folded-in"] 					= vertex2

			for neighbor in common_neighbors:
				neighbor_weight1 									= folded_graph.edges[vertex1, neighbor]["weight"]
				neighbor_weight2 									= folded_graph.edges[vertex2, neighbor]["weight"]
				neighbor_weight_new 								= neighbor_weight1 + neighbor_weight2
				folded_graph.edges[vertex2, neighbor]["weight"] 	= neighbor_weight_new

				if abs(neighbor_weight_new) < 1e-8:
					folded_graph.remove_edge(vertex2, neighbor)

				folded_graph.nodes[vertex2]["pos-weight"]			+= max(neighbor_weight_new, 0) - max(neighbor_weight2, 0)
				folded_graph.nodes[vertex2]["neg-weight"]			+= min(neighbor_weight_new, 0) - min(neighbor_weight2, 0)


			neighbors_of_vertex1 			= set(folded_graph.neighbors(vertex1))
			remaining_neighbors_vertex1		= list(neighbors_of_vertex1 - set(common_neighbors + [vertex2]) )

			for neighbor in remaining_neighbors_vertex1:
				folded_graph.add_edge(neighbor, vertex2)
				neighbor_weight 									= folded_graph.edges[vertex1, neighbor]["weight"]
				folded_graph.edges[vertex2, neighbor]["weight"] 	= neighbor_weight

				folded_graph.nodes[vertex2]["pos-weight"]			+= max(neighbor_weight, 0)
				folded_graph.nodes[vertex2]["neg-weight"]			+= min(neighbor_weight, 0)


			if self.fixed_vertex == vertex1:
				self.fixed_vertex 									= vertex2
				self.graph.nodes[vertex2]["partition"] 				= self.graph.nodes[vertex1]["partition"] 
				folded_graph.nodes[vertex2]["partition"] 			= self.graph.nodes[vertex1]["partition"]

			folded_graph.remove_node(vertex1)

	self.preprocessing_end_time				= time()
	self.preprocessing_total_time 			= self.preprocessing_end_time  - self.preprocessing_start_time


	return (is_folded, folded_graph)




#================================================================================================
# Update the parent of the current node in the decomposition tree with the folding operation applied
#================================================================================================
def update_parent_fold(self):

	self.parent.preprocessing_start_time			= time()

	for vertex in self.parent.vertices:
		folded_vertex								= vertex
		folded_vertex_temp 							= self.parent.graph.nodes[vertex]["folded-in"]

		while(folded_vertex_temp != -1):
			folded_vertex 							= folded_vertex_temp
			folded_vertex_temp 						= self.parent.graph.nodes[folded_vertex_temp]["folded-in"]

		self.parent.graph.nodes[vertex]["partition"]		= self.graph.nodes[folded_vertex]["partition"]

	self.parent.upper_bound 						= self.upper_bound
	
	self.parent.preprocessing_end_time				= time()
	self.parent.preprocessing_total_time 			+= self.parent.preprocessing_end_time  - self.parent.preprocessing_start_time

	


