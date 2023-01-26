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


		neighbors_of_vertex1 		= set(folded_graph.neighbors(vertex1) )
		neighbors_of_vertex2 		= set(folded_graph.neighbors(vertex2) )

		common_neighbors 			= list( neighbors_of_vertex1.intersection(neighbors_of_vertex2) )
		

		# min_common_weights 			= {neighbor: min(folded_graph.edges[vertex1, neighbor]["weight"], folded_graph.edges[vertex2, neighbor]["weight"]) for neighbor in common_neighbors}
		
		common_neighbors_weight		= {neighbor: (folded_graph.edges[vertex1, neighbor]["weight"], folded_graph.edges[vertex2, neighbor]["weight"]) for neighbor in common_neighbors}

		min_common_weights 			= {neighbor: (weight1 if abs(weight1) <= abs(weight2) else weight2) for (neighbor, (weight1, weight2) ) in common_neighbors_weight.items()}
		beta1_common_weights		= sum( (abs(weight1 - min_common_weights[neighbor]) if np.sign(weight1 * weight2) == -1 else 0) for (neighbor, (weight1, weight2) ) in common_neighbors_weight.items())
		beta2_common_weights		= sum( (abs(weight2 - min_common_weights[neighbor]) if np.sign(weight1 * weight2) == -1 else 0) for (neighbor, (weight1, weight2) ) in common_neighbors_weight.items())


		#----------------------------------------------------------------------------------------
		# vertices should satisfy the necessary and sufficient conditions
		#----------------------------------------------------------------------------------------	
		total_common_pos_weights	= sum(weight for weight in min_common_weights.values() if weight > 0)
		total_common_neg_weights	= sum(weight for weight in min_common_weights.values() if weight < 0)
		max_pos_diff_neg_weights	= max(folded_graph.nodes[vertex1]["pos-weight"] - folded_graph.nodes[vertex1]["neg-weight"] + beta1_common_weights, folded_graph.nodes[vertex2]["pos-weight"] - folded_graph.nodes[vertex2]["neg-weight"] + beta2_common_weights) 
		

		#----------------------------------------------------------------------------------------
		# Check the necessary condition
		#----------------------------------------------------------------------------------------
		neg_weight_of_edge 			= ( (2 if self.num_partitions == 2 else 1.5) * (min(folded_graph.edges[vertex1, vertex2]["weight"], 0) ) ) if folded_graph.has_edge(vertex1, vertex2) else 0
		can_be_folded 				= (total_common_pos_weights - total_common_neg_weights - neg_weight_of_edge >= 0.5 * max_pos_diff_neg_weights)
		
		can_be_folded_twin 			= (total_common_pos_weights - total_common_neg_weights - neg_weight_of_edge >= max_pos_diff_neg_weights)


		#----------------------------------------------------------------------------------------
		# Check the sufficient condition
		#----------------------------------------------------------------------------------------
		if can_be_folded == True and can_be_folded_twin == False:
	
			with Env(empty=True) as env:
				env.setParam('LogToConsole', 0)
				# sys.stdout.write(2*"\033[F\033[K")
				env.start()
				with Model(env=env) as model:

					t 					= model.addVars(common_neighbors, self.partitions, vtype=GRB.BINARY)
					q					= model.addVars(self.partitions, lb=- GRB.INFINITY, vtype=GRB.CONTINUOUS)
					
					model.setObjective(q[1] - q[0], GRB.MINIMIZE)

					for partition in self.partitions:
						model.addConstr(q[partition] == quicksum(min_common_weight * t[vertex, partition] for (vertex, min_common_weight) in min_common_weights.items()) )

					for partition in self.partitions[:-1]:
						model.addConstr(q[partition] <= q[partition + 1])

					for vertex in common_neighbors:
						model.addConstr( quicksum(t[vertex, partition] for partition in self.partitions) == 1)

					model.Params.Threads		= 1
					model.Params.LogToConsole 	= 0
					model.update()
					model.optimize()

					alpha_star					= model.objVal

					
					can_be_folded 				= (total_common_pos_weights - total_common_neg_weights - neg_weight_of_edge >= max_pos_diff_neg_weights - alpha_star)
	
		#-----------------------------------------------------------------------------------------
		# Fold vertices (graph.nodes[v1][folded-in] = v2  =>  v1 is removed and v2 is updated)
		#-----------------------------------------------------------------------------------------
		if can_be_folded == True:
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

	


