import networkx as nx 
import itertools
from copy import deepcopy


#================================================================================================
# Convert a binary string obtained from the quantum circuit to a binary solution 
#================================================================================================
def convert_string_sol_to_sorted_sol(self, string_solution):
	solution 			= [int(num) for num in string_solution]
	sorted_solution 	= {vertex:[0 for partition in self.partitions] for vertex in self.vertices} 

	sorted_solution[self.fixed_vertex] 		= [(1 if self.graph.nodes[self.fixed_vertex]["partition"] == partition else 0) for partition in self.partitions]

	if self.Params.Method == "QUBO" or self.Params.Method == "PUBO":
		for (v_ind, vertex) in enumerate(self.reduced_vertices):
			sorted_solution[vertex] 	= [solution[v_ind * self.num_partitions + partition] for partition in self.partitions]
	
	elif self.Params.Method == "R-QUBO":
		for (v_ind, vertex) in enumerate(self.reduced_vertices):
			sorted_solution[vertex] 	= [solution[v_ind * (self.num_partitions - 1) + partition] if partition != self.partitions[-1] else 0 for partition in self.partitions]
			sorted_solution[vertex][-1] = 1 - sum(sorted_solution[vertex]) 

	return sorted_solution



#================================================================================================
# Improve the resulting binary solution of QAOA and make it feasible in case of infeasibility 
#================================================================================================
def make_sol_feasible(self, solution):
	for vertex in self.reduced_vertices:
		num_assigned_partition 					= sum(solution[vertex])
		
		if num_assigned_partition != 1:
			allowed_partitions 					= [partition for partition in self.partitions if solution[vertex][partition] == 1] if num_assigned_partition != 0 else deepcopy(self.partitions)
			weight_partition_based_neighbors 	= {partition: sum([self.graph.edges[vertex, neighbor]["weight"] for neighbor in self.graph.neighbors(vertex) if solution[neighbor][partition] == 1]) for partition in allowed_partitions}
			selected_partition 					= min(weight_partition_based_neighbors, key=weight_partition_based_neighbors.get)
			solution[vertex] 					= [1 if partition == selected_partition else 0 for partition in self.partitions]


	return solution



#================================================================================================
# Calculate the objective function from the binary solution obtained from the quantum circuit
#================================================================================================
def cal_obj_from_sol(self, sorted_solution):
		
	#--------------------------------------------------------------------------------------------
	# Calculate the objective value of the sorted solution (in BQO format) 
	#--------------------------------------------------------------------------------------------
	objective_value 			= sum(self.graph.edges[edge[0], edge[1]]["weight"] for edge in self.graph.edges) 
	
	for edge in self.edges:
		objective_value 		-= self.graph.edges[edge[0], edge[1]]["weight"] * \
									sum([sorted_solution[edge[0]][partition] * sorted_solution[edge[1]][partition] for partition in self.partitions])

	#--------------------------------------------------------------------------------------------
	# Calculate the penalty term in the unconstrained binary optimization function
	#--------------------------------------------------------------------------------------------
	for vertex in self.vertices:
		penalty_term 			= 0
		
		if self.Params.Method == "QUBO":
			penalty_term 		= sum(sorted_solution[vertex]) - 1
			penalty_term 		= penalty_term * penalty_term

		elif self.Params.Method == "PUBO":
			penalty_term 		= 1
			for partition in self.partitions:
				penalty_term 	= penalty_term * (1 - sorted_solution[vertex][partition]) 

		elif self.Params.Method == "R-QUBO":
			for (partition1, partition2) in itertools.combinations(self.partitions[:-1], 2):
				penalty_term 	= penalty_term + sorted_solution[vertex][partition1] * sorted_solution[vertex][partition2]

		objective_value 		-= self.penalty_coef[vertex] * penalty_term

	return objective_value
