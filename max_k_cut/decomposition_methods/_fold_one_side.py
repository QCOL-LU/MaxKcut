import networkx as nx 
from gurobipy import *
import itertools
from copy import deepcopy
from time import time
import numpy as np

def folded_subgraph_solver(self, graph, folded_vertices, add_const=False):
	pass



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

		if vertex_pair_ditance[(vertex1, vertex2)] > 2: 
			continue


		if (not folded_graph.has_node(vertex1)) or (not folded_graph.has_node(vertex2)): continue


		neighbors_of_vertex1 		= set(folded_graph.neighbors(vertex1) ) - set([vertex2])
		neighbors_of_vertex2 		= set(folded_graph.neighbors(vertex2) ) - set([vertex1])

		
		beta 	= folded_graph.edges[vertex1, vertex2]["weight"] if folded_graph.has_edge(vertex1, vertex2) else 0

		if neighbors_of_vertex2 <= neighbors_of_vertex1:		# N(v1) is a subset of N(v2) 
			vertex1, vertex2 		= vertex2, vertex1
			neighbors_of_vertex1, neighbors_of_vertex2 	= neighbors_of_vertex2, neighbors_of_vertex1

		elif not (neighbors_of_vertex1 <= neighbors_of_vertex2): continue

				
		vertex1_neighbors_weight	= {neighbor: folded_graph.edges[vertex1, neighbor]["weight"] for neighbor in neighbors_of_vertex1}

		beta 						+= sum(abs(folded_graph.edges[vertex2, neighbor]["weight"]) for neighbor in (neighbors_of_vertex2 - neighbors_of_vertex1)) \
										+ sum(abs(folded_graph.edges[vertex2, neighbor]["weight"] - folded_graph.edges[vertex1, neighbor]["weight"]) for neighbor in neighbors_of_vertex1)

		#----------------------------------------------------------------------------------------
		# Check the necessary condition
		#----------------------------------------------------------------------------------------
		if (self.num_partitions > 2) and (beta > min(vertex1_neighbors_weight.values(),key=abs)): continue


		
		#----------------------------------------------------------------------------------------
		# Check the sufficient condition
		#----------------------------------------------------------------------------------------	
		with Env(empty=True) as env:
			env.setParam('LogToConsole', 0)
			# sys.stdout.write(2*"\033[F\033[K")
			env.start()
			with Model(env=env) as model:

				t 					= model.addVars(vertex1_neighbors_weight.keys(), self.partitions, vtype=GRB.BINARY)
				q					= model.addVars(self.partitions, lb=- GRB.INFINITY, vtype=GRB.CONTINUOUS)
				z 					= model.addVars(self.partitions, vtype=GRB.BINARY)
				
				my_list 			= sorted([abs(val) for val in vertex1_neighbors_weight.values()])
				temp 				= [(elm % 1)* 0.1 for elm in my_list]
				small_m 			= min(temp)
				big_M 				= sum(temp)
				
				model.addConstrs(q[partition] == quicksum(weight * t[neighbor, partition] 
								for (neighbor, weight) in vertex1_neighbors_weight.items()) for partition in self.partitions)
				
				model.addConstrs(q[partition] <= q[partition + 1] for partition in self.partitions[:-1])
				
				model.addConstrs(quicksum(t[vertex, partition] for partition in self.partitions) == 1 for vertex in vertex1_neighbors_weight)

				model.addConstr(quicksum(z[partition] for partition in self.partitions) == 1)

				for partition in self.partitions[1:]:
					model.addConstr( small_m * quicksum(z[inner_partition] for inner_partition in self.partitions[:partition + 1]) <= q[partition] - q[0])
					model.addConstr( big_M * quicksum(z[inner_partition] for inner_partition in self.partitions[:partition + 1]) >= q[partition] - q[0])
					model.addConstr( q[partition] - q[0] <= beta - small_m + big_M * (1 - z[partition]) ) 

				model.Params.Threads		= 1
				model.Params.LogToConsole 	= 0
				model.update()
				model.optimize()
				
				can_be_folded 				= not (model.status == GRB.Status.INFEASIBLE)
	
		#-----------------------------------------------------------------------------------------
		# Fold vertices (graph.nodes[v1][folded-in] = v2  =>  v1 is removed and v2 is updated)
		#-----------------------------------------------------------------------------------------
		if can_be_folded == True:
			is_folded 						= True
			self.applied_operation  		= "fold"

			self.graph.nodes[vertex1]["folded-in"] 					= vertex2

			for neighbor in vertex1_neighbors_weight:
				neighbor_weight1 									= folded_graph.edges[vertex1, neighbor]["weight"]
				neighbor_weight2 									= folded_graph.edges[vertex2, neighbor]["weight"]
				neighbor_weight_new 								= neighbor_weight1 + neighbor_weight2
				folded_graph.edges[vertex2, neighbor]["weight"] 	= neighbor_weight_new

				if abs(neighbor_weight_new) < 1e-8:
					folded_graph.remove_edge(vertex2, neighbor)

				folded_graph.nodes[vertex2]["pos-weight"]			+= max(neighbor_weight_new, 0) - max(neighbor_weight2, 0)
				folded_graph.nodes[vertex2]["neg-weight"]			+= min(neighbor_weight_new, 0) - min(neighbor_weight2, 0)



			# remaining_neighbors_vertex1		= set(folded_graph.neighbors(vertex2) ) - neighbors_of_vertex1

			# for neighbor in remaining_neighbors_vertex1:
			# 	folded_graph.add_edge(neighbor, vertex2)
			# 	neighbor_weight 									= folded_graph.edges[vertex1, neighbor]["weight"]
			# 	folded_graph.edges[vertex2, neighbor]["weight"] 	= neighbor_weight

			# 	folded_graph.nodes[vertex2]["pos-weight"]			+= max(neighbor_weight, 0)
			# 	folded_graph.nodes[vertex2]["neg-weight"]			+= min(neighbor_weight, 0)


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

	


