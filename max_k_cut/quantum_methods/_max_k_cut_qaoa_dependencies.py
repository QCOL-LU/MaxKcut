import networkx as nx 
import itertools
from copy import deepcopy
import sys 
import numpy as np
from qiskit.circuit.library.standard_gates import *



#================================================================================================
# Convert a binary string obtained from the quantum circuit to a binary solution 
#================================================================================================
def convert_string_sol_to_sorted_sol(self, string_solution):
	solution 			= [int(num) for num in string_solution]
	sorted_solution 	= {vertex:[0 for partition in self.partitions] for vertex in self.vertices} 

	sorted_solution[self.fixed_vertex] 		= [(1 if self.graph.nodes[self.fixed_vertex]["partition"] == partition else 0) for partition in self.partitions]

	rqubo_sorted_solution					= deepcopy(sorted_solution)

	if self.Params.Method == "QUBO" or self.Params.Method == "PUBO":
		for (v_ind, vertex) in enumerate(self.reduced_vertices):
			sorted_solution[vertex] 	= [solution[v_ind * self.num_partitions + partition] for partition in self.partitions]
			rqubo_sorted_solution[vertex] 		= deepcopy(sorted_solution[vertex])
			rqubo_sorted_solution[vertex][-1] 	= 1 - sum(sorted_solution[vertex][:-1]) 
	
	elif self.Params.Method == "R-QUBO":
		for (v_ind, vertex) in enumerate(self.reduced_vertices):
			sorted_solution[vertex] 	= [solution[v_ind * (self.num_partitions - 1) + partition] if partition != self.partitions[-1] else 0 for partition in self.partitions]
			sorted_solution[vertex][-1] = 1 - sum(sorted_solution[vertex]) 

		rqubo_sorted_solution 			= deepcopy(sorted_solution)
	return sorted_solution, rqubo_sorted_solution



#================================================================================================
# Improve the resulting binary solution of QAOA and make it feasible in case of infeasibility 
#================================================================================================
def make_sol_feasible(self, solution):

	#--------------------------------------------------------------------------------------------
	# make a solution of the QUBO formulation feasible
	#--------------------------------------------------------------------------------------------
	if self.Params.Method == "QUBO":
		vertex_num_assigned_partitions 			= {vertex: sum(solution[vertex]) for vertex in self.reduced_vertices}
		I_1 									= [vertex for vertex, num_assigned_partitions in vertex_num_assigned_partitions.items() if num_assigned_partitions > 1]
		I_0 									= [vertex for vertex, num_assigned_partitions in vertex_num_assigned_partitions.items() if num_assigned_partitions == 0]
		
		#----------------------------------------------------------------------------------------
		# Loop over vertices with multiple assignments
		#----------------------------------------------------------------------------------------
		while bool(I_1):
			selected_vertex						= I_1[0]
			common_assignment 					= []

			for vertex in I_1:
				flag 							= True
				
				for partition in self.partitions:
					if solution[vertex][partition] != solution[selected_vertex][partition]:
						flag 					= False
						break

				if flag == True: common_assignment.append(vertex)


			edges_ell 							= [edge for edge in nx.edge_boundary(self.graph, common_assignment, self.vertices) if self.graph.edges[edge]["weight"] < 0]
			partitions_ell 						= [partition for partition in self.partitions if solution[selected_vertex][partition] == 1]
			
			weight_partition_based_neighbors 	= {partition: sum(self.graph.edges[edge]["weight"] for edge in edges_ell if solution[edge[0]][partition] == solution[edge[1]][partition]) for partition in partitions_ell}
			selected_partition 					= min(weight_partition_based_neighbors, key=weight_partition_based_neighbors.get)


			for vertex in common_assignment:
				for partition in partitions_ell:
					solution[vertex][partition] = 0 if partition != selected_partition else solution[vertex][partition]

				I_1.remove(vertex)
		

		#----------------------------------------------------------------------------------------
		# Loop over vertices with no assignment
		#----------------------------------------------------------------------------------------
		for vertex in I_0:
			weight_partition_based_neighbors 	= {partition: sum([self.graph.edges[vertex, neighbor]["weight"] for neighbor in self.graph.neighbors(vertex) if solution[neighbor][partition] == 1]) for partition in self.partitions}
			selected_partition 					= min(weight_partition_based_neighbors, key=weight_partition_based_neighbors.get)
			solution[vertex][selected_partition]= 1


	#--------------------------------------------------------------------------------------------
	# make a solution of the R-QUBO formulation feasible
	#--------------------------------------------------------------------------------------------
	elif self.Params.Method == "R-QUBO":
		all_vertices 			= {vertex: sum(solution[vertex][partition] for partition in self.partitions[:-1]) for vertex in self.reduced_vertices}
		M_ordering 				= {vertex: value for (vertex, value) in all_vertices.items() if value > 1}

		while  bool(M_ordering):
			selected_vertex 		= max(M_ordering, key=M_ordering.get)
			selected_value 			= M_ordering[selected_vertex]

			common_assignment 		= []

			for vertex in M_ordering.keys():
				flag 	= True
				
				for partition in self.partitions[:-1]:
					if solution[vertex][partition] != solution[selected_vertex][partition]:
						flag = False
						break

				if flag == True: common_assignment.append(vertex)


			edges_ell 				= nx.edge_boundary(self.graph, common_assignment, self.vertices)
			partitions_ell 			= [partition for partition in self.partitions[:-1] if solution[selected_vertex][partition] == 1] 
			

			weight_partition_based_neighbors 	= {partition: sum(self.graph.edges[edge]["weight"] for edge in edges_ell if solution[edge[0]][partition] == solution[edge[1]][partition]) for partition in partitions_ell}
			selected_partition 					= min(weight_partition_based_neighbors, key=weight_partition_based_neighbors.get)


			for vertex in common_assignment:
				for partition in partitions_ell:
					solution[vertex][partition] = 0 if partition != selected_partition else solution[vertex][partition]


				del M_ordering[vertex]
				solution[vertex][-1] = 1 - sum(solution[vertex][:-1])

	return solution



#================================================================================================
# Calculate the objective function from the binary solution obtained from the quantum circuit
#================================================================================================
def cal_obj_from_sol(self, sorted_solution):
		
	
	#--------------------------------------------------------------------------------------------
	# Calculate the objective value of the sorted solution (in BQO format) 
	#--------------------------------------------------------------------------------------------
	total_penalty 				= 0
	tight_total_penalty			= 0
	objective_value 			= sum(self.graph.edges[edge[0], edge[1]]["weight"] for edge in self.graph.edges) 
	
	for edge in self.edges:
		objective_value 		-= self.graph.edges[edge[0], edge[1]]["weight"] * \
									sum([sorted_solution[edge[0]][partition] * sorted_solution[edge[1]][partition] for partition in self.partitions])

	if self.Params.Method == "QUBO":
		vertex_num_assigned_partitions	= {vertex: sum(sorted_solution[vertex]) for vertex in self.reduced_vertices}
		constraint_violation			= sum(abs(vertex_num_assigned_partitions[vertex] - 1) for vertex in self.reduced_vertices) 

	elif self.Params.Method == "R-QUBO":
		vertex_num_assigned_partitions	= {vertex: sum(sorted_solution[vertex][partition] for partition in self.partitions[:-1]) for vertex in self.reduced_vertices}
		constraint_violation 			= sum(max(vertex_num_assigned_partitions[vertex] - 1, 0) for vertex in self.reduced_vertices)
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

		total_penalty 			+= self.penalty_coef[vertex] * penalty_term
		tight_total_penalty 	+= self.tight_penalty_coef[vertex] * penalty_term
	
	tight_objective_value		= objective_value - tight_total_penalty
	objective_value 			= objective_value - total_penalty
	

	return tight_objective_value, objective_value, total_penalty, constraint_violation


#================================================================================================
# Calculate the average and the best objective values and the best solution 
#================================================================================================
def cal_avg_best_sol(self, angles):
	shots							= float(self.Params.QAOA_Num_Shots)

	sum_obj_value					= 0.0
	sum_tight_obj_value				= 0.0
	sum_rqubo_tight_obj_value		= 0.0

	sum_total_penalty 				= 0.0
	sum_rqubo_total_penalty 		= 0.0
	sum_constraint_violation 		= 0.0

	sum_pure_feasible_obj_value 	= 0.0
	pure_feasible_count 			= 0.0



	self.histogram 					= {num: 0 for num in range(len(self.graph.edges))}

	for (string_solution, count) in self.counts.items():
		sorted_solution, rqubo_sorted_solution			= convert_string_sol_to_sorted_sol(self, string_solution)
		tight_obj_value, obj_value, total_penalty, constraint_violation 	= cal_obj_from_sol(self, sorted_solution)
		temp_method 				= self.Params.Method
		self.Params.Method 			= "R-QUBO"

		rqubo_tight_obj_value,  rqubo_tight_total_penalty		= (tight_obj_value, total_penalty )if temp_method == "R-QUBO" else cal_obj_from_sol(self, rqubo_sorted_solution)[1:3]
		self.Params.Method 			= temp_method


		count 					 	= int(count * shots) if count < 1 else count

		if total_penalty == 0.0:
			pure_feasible_count 		= pure_feasible_count + count
			sum_pure_feasible_obj_value = sum_pure_feasible_obj_value + obj_value * count


		
		sum_obj_value				= sum_obj_value + count * obj_value
		sum_tight_obj_value 		= sum_tight_obj_value + count * tight_obj_value
		sum_rqubo_tight_obj_value 	= sum_rqubo_tight_obj_value + count * rqubo_tight_obj_value

		sum_total_penalty 			= sum_total_penalty + count * total_penalty
		sum_rqubo_total_penalty 	= sum_rqubo_total_penalty + count *rqubo_tight_total_penalty
		sum_constraint_violation 	= sum_constraint_violation + count * constraint_violation

		self.histogram[obj_value]	= self.histogram.get(obj_value, 0) + 100 * count/shots

		if self.qaoa_best_obj_value < obj_value:
			self.qaoa_best_obj_value 	= obj_value
			self.qaoa_best_solution 	= deepcopy(sorted_solution)


	

	self.qaoa_avg_obj_value 				= sum_obj_value / shots
	self.qaoa_avg_tight_obj_value			= sum_tight_obj_value / shots
	self.qaoa_avg_rqubo_tight_obj_value		= sum_rqubo_tight_obj_value / shots

	self.qaoa_avg_total_penalty  			= sum_total_penalty / shots
	self.qaoa_avg_rqubo_total_penalty  		= sum_rqubo_total_penalty / shots
	self.qaoa_avg_constraint_violation  	= sum_constraint_violation / shots

	self.qaoa_avg_pure_feasible_obj_value 	= sum_pure_feasible_obj_value / pure_feasible_count if pure_feasible_count != 0 else "NA"
	self.qaoa_pure_feasible_percentage		= pure_feasible_count / shots * 100


	biased_sample 							= {value - self.qaoa_avg_obj_value: count for (value, count) in self.histogram.items()}
	biased_sample_2 						= sum(count * value**2 for (value, count) in  biased_sample.items() ) / shots
	biased_sample_3 						= sum(count * value**3 for (value, count) in  biased_sample.items() ) / shots

	self.qaoa_skewness 						= np.sqrt(shots * (shots - 1)) / (shots - 1) * (biased_sample_3) / np.sqrt(biased_sample_2)**3

	if self.qaoa_best_avg_obj_value < self.qaoa_avg_obj_value:
		self.qaoa_best_avg_obj_value 				= self.qaoa_avg_obj_value
		self.qaoa_best_avg_tight_obj_value			= self.qaoa_avg_tight_obj_value
		self.qaoa_best_avg_rqubo_tight_obj_value	= self.qaoa_avg_rqubo_tight_obj_value

		self.qaoa_best_avg_total_penalty 			= self.qaoa_avg_total_penalty
		self.qaoa_best_avg_rqubo_total_penalty 		= self.qaoa_avg_rqubo_total_penalty
		self.qaoa_best_avg_constraint_violation 	= self.qaoa_avg_constraint_violation

		self.qaoa_best_avg_pure_feasible_obj_value	= self.qaoa_avg_pure_feasible_obj_value
		self.qaoa_best_pure_feasible_percentage 	= self.qaoa_pure_feasible_percentage

		self.qaoa_best_skewness 					= self.qaoa_skewness
		self.best_angles							= angles
		self.best_histogram 						= self.histogram

		self.best_counts 							= self.counts

	else:
		self.qaoa_best_ev_improved 					= False



		
	

	self.qaoa_std_obj_value 					= np.sqrt(sum([(key - self.qaoa_avg_obj_value)**2 * value for (key, value) in self.histogram.items()])/(shots - 1))	



#================================================================================================
# Calculate the average and the best objective values and the best solution 
#================================================================================================
def cal_avg_best_sol_feasible(self):
	shots							= float(self.Params.QAOA_Num_Shots)

	feasible_sum_obj_value			= 0.0
	self.feasible_histogram 		= {num: 0 for num in range(len(self.graph.edges))}


	for (string_solution, count) in self.best_counts.items():
		sorted_solution, rqubo_sorted_solution 	= convert_string_sol_to_sorted_sol(self, string_solution)

		feasible_solution						= make_sol_feasible(self, sorted_solution)
		tight_obj_value, feasible_obj_value,	total_penalty, constraint_violation		= cal_obj_from_sol(self, feasible_solution)
		feasible_sum_obj_value					= feasible_sum_obj_value + count * feasible_obj_value


		self.feasible_histogram[feasible_obj_value]	= self.feasible_histogram.get(feasible_obj_value, 0) + 100 * count/shots


		if self.qaoa_feasible_best_obj_value < feasible_obj_value:
			self.qaoa_feasible_best_obj_value 	= feasible_obj_value
			self.qaoa_feasible_best_solution 	= deepcopy(feasible_solution)


	self.qaoa_feasible_avg_obj_value 		= feasible_sum_obj_value / shots



	if self.qaoa_feasible_best_avg_obj_value < self.qaoa_feasible_avg_obj_value:
		self.qaoa_feasible_best_avg_obj_value	= self.qaoa_feasible_avg_obj_value
		self.feasible_best_histogram 			= self.feasible_histogram

	
	
	self.qaoa_feasible_std_obj_value 			= np.sqrt(sum([(key - self.qaoa_feasible_avg_obj_value)**2 * value for (key, value) in self.feasible_histogram.items()])/(shots - 1))
	




#================================================================================================
# Apply gate exp(i * angle * (II - ZZ)) on the QAOA circuit
#================================================================================================
def gate_i_zz(self, angle, qubits):
	[q1, q2] 	= qubits
	U1 			= U1Gate(angle)
	CU1 		= U1Gate(-2*angle).control(1)
	self.qaoa_circuit.append(CU1, [q1, q2])
	self.qaoa_circuit.append(U1, [q1])
	self.qaoa_circuit.append(U1, [q2])

#================================================================================================
# Apply gate exp(i * angle * (I - Z)) on the QAOA circuit
#================================================================================================
def gate_i_z_1(self, angle, qubit):
	U1 		= U1Gate(angle)
	self.qaoa_circuit.append(U1, qubit)

#================================================================================================
# Apply gate exp(i * angle * (I - Z)(I - Z)) on the QAOA circuit
#================================================================================================
def gate_i_z_2(self, angle, qubits):
	CU1 		= U1Gate(angle).control(1)
	self.qaoa_circuit.append(CU1, qubits)




