import networkx as nx 
import sys 
import itertools
import time
import numpy as np


from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, execute, Aer
from qiskit.extensions import HamiltonianGate

from quantum_methods._max_k_cut_qaoa_dependencies import cal_obj_from_sol, convert_string_sol_to_sorted_sol, make_sol_feasible


#================================================================================================
# Calculate the objective function from the binary solution obtained from quantum circuit
#================================================================================================
def qaoa_expected_value(self, angles):

	#--------------------------------------------------------------------------------------------
	# Parameters setting 
	#--------------------------------------------------------------------------------------------
	
	if self.Params.Method == "QUBO":
		self.penalty_coef	= {vertex:  self.Params.Penalty_Increase * max(self.graph.nodes[vertex]["pos-weight"]/self.num_partitions, - self.graph.nodes[vertex]["neg-weight"] / 2) for vertex in self.vertices}
		num_qubits 			= (self.num_vertices - 1) * self.num_partitions
		qubit_map 			= {vertex:[v_ind * self.num_partitions + partition for partition in self.partitions] for (v_ind, vertex) in enumerate(self.reduced_vertices)}
	
	elif self.Params.Method == "PUBO":
		self.penalty_coef	= {vertex: self.Params.Penalty_Increase * max(self.graph.nodes[vertex]["pos-weight"]/self.num_partitions, -(self.num_partitions - 1) * self.graph.nodes[vertex]["neg-weight"]) for vertex in self.vertices}
		num_qubits 			= (self.num_vertices - 1) * self.num_partitions
		qubit_map 			= {vertex:[v_ind * self.num_partitions + partition for partition in self.partitions] for (v_ind, vertex) in enumerate(self.reduced_vertices)}

	elif self.Params.Method == "R-QUBO":
		self.penalty_coef	= {vertex: self.Params.Penalty_Increase * (self.graph.nodes[vertex]["pos-weight"] -  self.graph.nodes[vertex]["neg-weight"])  for vertex in self.vertices}
		num_qubits 			= (self.num_vertices - 1) * (self.num_partitions - 1)
		qubit_map 			= {vertex:[v_ind * (self.num_partitions - 1) + partition for partition in self.partitions[:-1]] for (v_ind, vertex) in enumerate(self.reduced_vertices)}
	
	all_qubits 				= [qubit for qubit in range(num_qubits)]


	#--------------------------------------------------------------------------------------------
	# Quantum circuit initialization (QUBO, PUBO, and R-QUBO)
	#--------------------------------------------------------------------------------------------
	self.qaoa_circuit 		= QuantumCircuit(num_qubits, num_qubits)
	self.qaoa_circuit.h(all_qubits)


	for level in range(self.Params.QAOA_Num_Levels):
		gamma 			= angles[level]
		beta 			= angles[self.Params.QAOA_Num_Levels + level]

		#---------------------------------------------------------------------------------------
		# Phase separation operator - QUBO & PUBO (objective function of the BQO)
		#---------------------------------------------------------------------------------------
		if self.Params.Method == "QUBO" or self.Params.Method == "PUBO":

			#-----------------------------------------------------------------------------------
			# Unitary operators associated to the objective of BQO and non-incident edges of the fixed vertex (QUBO & PUBO)
			#-----------------------------------------------------------------------------------
			for edge in self.edges:
				if edge[0] != self.fixed_vertex and edge[1] != self.fixed_vertex:
					hamiltonian_matrix 		= np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])
					unitary_operator 		= HamiltonianGate(hamiltonian_matrix, - gamma * self.graph.edges[edge[0], edge[1]]["weight"])
					
					for partition in self.partitions:
						qubits 				= [qubit_map[edge[0]][partition], qubit_map[edge[1]][partition]]
						self.qaoa_circuit.append(unitary_operator, qubits)

			#-----------------------------------------------------------------------------------
			# Unitary operators associated to the objective of BQO and incident edges of the fixed vertex (QUBO & PUBO)
			#-----------------------------------------------------------------------------------
			partition 				= self.graph.nodes[self.fixed_vertex]["partition"]
			hamiltonian_matrix 		= np.array([[0,0],[0,1]])

			for neighbor in self.graph.neighbors(self.fixed_vertex):
				edge_weight 		= self.graph.edges[self.fixed_vertex, neighbor]["weight"]
				unitary_operator 	= HamiltonianGate(hamiltonian_matrix, - gamma * edge_weight)
				qubit 				= [qubit_map[neighbor][partition]]
				self.qaoa_circuit.append(unitary_operator, qubit)


		#---------------------------------------------------------------------------------------
		# Phase separation operator - QUBO (penalty terms of constraint violation)
		#---------------------------------------------------------------------------------------
		if self.Params.Method == "QUBO":
			for vertex in self.reduced_vertices:

				#-------------------------------------------------------------------------------
				# Unitary operators associated to pair of partitions (QUBO)
				#-------------------------------------------------------------------------------
				hamiltonian_matrix 	= np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])
				unitary_operator 	= HamiltonianGate(hamiltonian_matrix, - 2 * gamma * self.penalty_coef[vertex])

				for (partition1, partition2) in itertools.combinations(self.partitions, 2):
					qubits 			= [qubit_map[vertex][partition1], qubit_map[vertex][partition2]]
					self.qaoa_circuit.append(unitary_operator, qubits)

				#-------------------------------------------------------------------------------
				# Unitary operators associated to each of partitions (QUBO)
				#-------------------------------------------------------------------------------
				hamiltonian_matrix 	= np.array([[0,0],[0,1]])
				unitary_operator 	= HamiltonianGate(hamiltonian_matrix, gamma * self.penalty_coef[vertex])

				for partition in self.partitions:
					qubit 			= [qubit_map[vertex][partition]]
					self.qaoa_circuit.append(unitary_operator, qubit)


		#---------------------------------------------------------------------------------------
		# Phase separation operator - PUBO (penalty terms of constraint violation)
		#---------------------------------------------------------------------------------------
		if self.Params.Method == "PUBO":
			hamiltonian_matrix 		= np.zeros((2**self.num_partitions, 2**self.num_partitions))
			hamiltonian_matrix[0,0]	= 1
	
			for vertex in self.reduced_vertices:
				unitary_operator 	= HamiltonianGate(hamiltonian_matrix, - gamma * self.penalty_coef[vertex])
				qubits 				= [qubit_map[vertex][partition] for partition in self.partitions]
				self.qaoa_circuit.append(unitary_operator, qubits)


		#---------------------------------------------------------------------------------------
		# Phase separation operator - R-QUBO 
		#---------------------------------------------------------------------------------------
		if self.Params.Method == "R-QUBO":

			#-----------------------------------------------------------------------------------
			# Unitary operators associated to the objective of R-BQO and non-incident edges of the fixed vertex (R-QUBO)
			#-----------------------------------------------------------------------------------
			for edge in self.edges:
				if edge[0] != self.fixed_vertex and edge[1] != self.fixed_vertex:

					#---------------------------------------------------------------------------
					# Unitary operators associated to edge on xor operation (R-QUBO)
					#---------------------------------------------------------------------------
					hamiltonian_matrix 		= np.array([[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,0]])
					edge_weight 			= self.graph.edges[edge[0], edge[1]]["weight"]
					unitary_operator 		= HamiltonianGate(hamiltonian_matrix, gamma * edge_weight)
					
					for partition in self.partitions[:-1]:
						qubits 				= [qubit_map[edge[0]][partition], qubit_map[edge[1]][partition]]
						self.qaoa_circuit.append(unitary_operator, qubits)


					#---------------------------------------------------------------------------
					# Unitary operators associated to pair of partitions (R-QUBO)
					#---------------------------------------------------------------------------
					hamiltonian_matrix 		= np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])
					unitary_operator 		= HamiltonianGate(hamiltonian_matrix, - gamma * edge_weight)

					for (partition1, partition2) in itertools.combinations(self.partitions[:-1], 2):
						qubits1 			= [qubit_map[edge[0]][partition1], qubit_map[edge[1]][partition2]]
						self.qaoa_circuit.append(unitary_operator, qubits1)

						qubits2 			= [qubit_map[edge[1]][partition1], qubit_map[edge[0]][partition2]]
						self.qaoa_circuit.append(unitary_operator, qubits2)


			#-----------------------------------------------------------------------------------
			# Unitary operators associated to the objective of R-BQO and incident edges of the fixed vertex (R-QUBO)
			#-----------------------------------------------------------------------------------
			fixed_vertex_partition 		= self.graph.nodes[self.fixed_vertex]["partition"]
			hamiltonian_matrix 			= np.array([[1,0],[0,0]])

			for neighbor in self.graph.neighbors(self.fixed_vertex):

				edge_weight 			= self.graph.edges[neighbor, self.fixed_vertex]["weight"]
				
				
				if fixed_vertex_partition != self.partitions[-1]:
					unitary_operator	= HamiltonianGate(hamiltonian_matrix, - gamma * edge_weight)
					qubit 				= [qubit_map[neighbor][fixed_vertex_partition]]
					self.qaoa_circuit.append(unitary_operator, qubit)

				else:
					unitary_operator	= HamiltonianGate(hamiltonian_matrix, gamma * edge_weight)
					for partition in self.partitions[:-1]:
						qubit 			= [qubit_map[neighbor][partition]]
						self.qaoa_circuit.append(unitary_operator, qubit)



			#-----------------------------------------------------------------------------------
			# Unitary operators associated to penalty terms of constraint violation (R-QUBO)
			#-----------------------------------------------------------------------------------
			for vertex in self.reduced_vertices:
				hamiltonian_matrix 		= np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])
				unitary_operator 		= HamiltonianGate(hamiltonian_matrix, - gamma * self.penalty_coef[vertex])

				for (partition1, partition2) in itertools.combinations(self.partitions[:-1], 2):
					qubits 				= [qubit_map[vertex][partition1], qubit_map[vertex][partition2]]
					self.qaoa_circuit.append(unitary_operator, qubits)


		#---------------------------------------------------------------------------------------
		# Mixing operator (QUBO, PUBO, and R-QUBO)
		#---------------------------------------------------------------------------------------
		self.qaoa_circuit.barrier()
		self.qaoa_circuit.rx(2 * beta, all_qubits)
		self.qaoa_circuit.barrier()


	#--------------------------------------------------------------------------------------------
	# Measurement (QUBO, PUBO, and R-QUBO)
	#--------------------------------------------------------------------------------------------
	self.qaoa_circuit.measure(all_qubits, all_qubits)


	backend 			= Aer.get_backend("qasm_simulator") #qasm_simulator
	shots				= self.Params.QAOA_Num_Shots

	
	simulate 			= execute(self.qaoa_circuit, backend=backend, shots=shots)
	qaoa_results		= simulate.result()
	
	self.counts 		= qaoa_results.get_counts()


	#--------------------------------------------------------------------------------------------
	# Calculate the average and the best objective values and the best solution 
	#--------------------------------------------------------------------------------------------
	best_obj_value 					= - sys.maxsize
	sum_obj_value					= 0.0
	self.histogram 					= {num: 0 for num in range(len(self.graph.edges))}

	feasible_best_obj_value 		= - sys.maxsize
	feasible_sum_obj_value			= 0.0
	self.feasible_histogram 		= {num: 0 for num in range(len(self.graph.edges))}


	for (string_solution, count) in self.counts.items():
		sorted_solution 			= convert_string_sol_to_sorted_sol(self, string_solution)
		obj_value 					= cal_obj_from_sol(self, sorted_solution)
		sum_obj_value				= sum_obj_value + count * obj_value

		self.histogram[obj_value]	= self.histogram.get(obj_value, 0) + 100 * count/shots


		feasible_solution			= make_sol_feasible(self, sorted_solution)
		feasible_obj_value			= cal_obj_from_sol(self, feasible_solution)
		feasible_sum_obj_value		= feasible_sum_obj_value + count * feasible_obj_value

		self.feasible_histogram[obj_value]	= self.histogram.get(feasible_obj_value, 0) + 100 * count/shots


		if best_obj_value < obj_value:
			best_obj_value 			= obj_value
			best_solution			= sorted_solution

		if feasible_best_obj_value < feasible_obj_value:
			feasible_best_obj_value = feasible_obj_value
			feasible_best_solution	= feasible_solution


	if self.qaoa_best_obj_value < best_obj_value:
		self.qaoa_best_obj_value 	= best_obj_value
		self.qaoa_best_solution 	= best_solution

	if self.qaoa_feasible_best_obj_value < feasible_best_obj_value:
		self.qaoa_feasible_best_obj_value 	= feasible_best_obj_value
		self.qaoa_feasible_best_solution 	= feasible_best_solution


	self.qaoa_avg_obj_value 				= sum_obj_value / float(shots)

	self.qaoa_feasible_avg_obj_value 		= feasible_sum_obj_value / float(shots)


	if self.qaoa_best_avg_obj_value < self.qaoa_avg_obj_value:
		self.qaoa_best_avg_obj_value= self.qaoa_avg_obj_value
		self.best_angles			= angles
		self.best_histogram 		= self.histogram


	if self.qaoa_feasible_best_avg_obj_value < self.qaoa_feasible_avg_obj_value:
		self.qaoa_feasible_best_avg_obj_value	= self.qaoa_feasible_avg_obj_value
		self.feasible_best_angles				= angles
		self.feasible_best_histogram 			= self.feasible_histogram

	self.qaoa_circuit_depth						= self.qaoa_circuit.depth()
	

	self.qaoa_std_obj_value 					= np.sqrt(sum([(key - self.qaoa_avg_obj_value)**2 * value for (key, value) in self.histogram.items()])/(shots - 1))
	self.qaoa_feasible_std_obj_value 			= np.sqrt(sum([(key - self.qaoa_feasible_avg_obj_value)**2 * value for (key, value) in self.feasible_histogram.items()])/(shots - 1))
	

	#--------------------------------------------------------------------------------------------
	# Print the results of QAOA optimizer at each iteration  
	#--------------------------------------------------------------------------------------------
	if self.Params.QAOA_Optimize == True:

		self.print_qaoa_optimizer_iter()
		self.qaoa_iter 					= self.qaoa_iter + 1

	else:
		self.qoao_opt_end_time 			= time.time()
		self.qaoa_opt_total_time 		= self.qoao_opt_end_time - self.qoao_opt_start_time


	return - self.qaoa_avg_obj_value



