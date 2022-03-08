import networkx as nx 
import sys 
import itertools
import time
import numpy as np


from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, execute, Aer
from qiskit.extensions import HamiltonianGate
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import  depolarizing_error

from quantum_methods._max_k_cut_qaoa_dependencies import cal_obj_from_sol, convert_string_sol_to_sorted_sol, make_sol_feasible, cal_avg_best_sol, cal_avg_best_sol_feasible, gate_i_zz, gate_i_z_1, gate_i_z_2


#================================================================================================
# Calculate the objective function from the binary solution obtained from quantum circuit
#================================================================================================
def qaoa_expected_value(self, angles):

	#--------------------------------------------------------------------------------------------
	# Parameters setting 
	#--------------------------------------------------------------------------------------------
	if self.Params.Method == "QUBO":
		self.penalty_coef	= {vertex:  self.Params.Penalty_Increase * max(self.graph.nodes[vertex]["pos-weight"]/self.num_partitions, - self.graph.nodes[vertex]["neg-weight"] / 2) + 0.1 for vertex in self.vertices}
		num_qubits 			= (self.num_vertices - 1) * self.num_partitions
		qubit_map 			= {vertex:[v_ind * self.num_partitions + partition for partition in self.partitions] for (v_ind, vertex) in enumerate(self.reduced_vertices)}
	
		if self.Params.Naive == True:
			self.penalty_coef	= {vertex: (self.graph.nodes[vertex]["pos-weight"] -  self.graph.nodes[vertex]["neg-weight"]) for vertex in self.vertices}


	elif self.Params.Method == "PUBO":
		self.penalty_coef	= {vertex: self.Params.Penalty_Increase * max(self.graph.nodes[vertex]["pos-weight"]/self.num_partitions, -(self.num_partitions - 1) * self.graph.nodes[vertex]["neg-weight"]) for vertex in self.vertices}
		num_qubits 			= (self.num_vertices - 1) * self.num_partitions
		qubit_map 			= {vertex:[v_ind * self.num_partitions + partition for partition in self.partitions] for (v_ind, vertex) in enumerate(self.reduced_vertices)}

	elif self.Params.Method == "R-QUBO":
		self.penalty_coef	= {vertex: self.Params.Penalty_Increase * (self.graph.nodes[vertex]["pos-weight"] -  self.graph.nodes[vertex]["neg-weight"])  for vertex in self.vertices}
		num_qubits 			= (self.num_vertices - 1) * (self.num_partitions - 1)
		qubit_map 			= {vertex:[v_ind * (self.num_partitions - 1) + partition for partition in self.partitions[:-1]] for (v_ind, vertex) in enumerate(self.reduced_vertices)}
		
		if self.Params.Naive == True:
			self.penalty_coef	= {vertex: self.num_partitions * (self.graph.nodes[vertex]["pos-weight"] -  self.graph.nodes[vertex]["neg-weight"])  for vertex in self.vertices}


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
		# Phase separation operator - QUBO  (objective function of the BQO)
		#---------------------------------------------------------------------------------------
		if self.Params.Method == "QUBO":

			#-----------------------------------------------------------------------------------
			# Unitary operators associated to the objective of BQO and non-incident edges of the fixed vertex (QUBO)
			#-----------------------------------------------------------------------------------
			for edge in self.edges:
				if edge[0] != self.fixed_vertex and edge[1] != self.fixed_vertex:
										
					for partition in self.partitions:
						qubits 				= [qubit_map[edge[0]][partition], qubit_map[edge[1]][partition]]
						self.gate_i_z_2(gamma * self.graph.edges[edge[0], edge[1]]["weight"], qubits)

			#-----------------------------------------------------------------------------------
			# Unitary operators associated to the objective of BQO and incident edges of the fixed vertex (QUBO)
			#-----------------------------------------------------------------------------------
			partition 				= self.graph.nodes[self.fixed_vertex]["partition"]

			for neighbor in self.graph.neighbors(self.fixed_vertex):
				edge_weight 		= self.graph.edges[self.fixed_vertex, neighbor]["weight"]
				qubit 				= [qubit_map[neighbor][partition]]
				self.gate_i_z_1(gamma * edge_weight, qubit)


		#---------------------------------------------------------------------------------------
		# Phase separation operator - QUBO (penalty terms of constraint violation)
		#---------------------------------------------------------------------------------------
		if self.Params.Method == "QUBO":
			for vertex in self.reduced_vertices:

				#-------------------------------------------------------------------------------
				# Unitary operators associated to pair of partitions (QUBO)
				#-------------------------------------------------------------------------------
				for (partition1, partition2) in itertools.combinations(self.partitions, 2):
					qubits 			= [qubit_map[vertex][partition1], qubit_map[vertex][partition2]]
					self.gate_i_z_2(2 * gamma * self.penalty_coef[vertex], qubits)

				#-------------------------------------------------------------------------------
				# Unitary operators associated to each of partitions (QUBO)
				#-------------------------------------------------------------------------------
				for partition in self.partitions:
					qubit 			= [qubit_map[vertex][partition]]
					self.gate_i_z_1(- gamma * self.penalty_coef[vertex], qubit)





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
					edge_weight 			= self.graph.edges[edge[0], edge[1]]["weight"]					
					for partition in self.partitions[:-1]:
						qubits 				= [qubit_map[edge[0]][partition], qubit_map[edge[1]][partition]]
						self.gate_i_zz(- gamma * edge_weight, qubits)


					#---------------------------------------------------------------------------
					# Unitary operators associated to pair of partitions (R-QUBO)
					#---------------------------------------------------------------------------
					for (partition1, partition2) in itertools.combinations(self.partitions[:-1], 2):
						qubits1 			= [qubit_map[edge[0]][partition1], qubit_map[edge[1]][partition2]]
						qubits2 			= [qubit_map[edge[1]][partition1], qubit_map[edge[0]][partition2]]
						
						self.gate_i_z_2(gamma * edge_weight, qubits1)
						self.gate_i_z_2(gamma * edge_weight, qubits2)


			#-----------------------------------------------------------------------------------
			# Unitary operators associated to the objective of R-BQO and incident edges of the fixed vertex (R-QUBO)
			#-----------------------------------------------------------------------------------
			fixed_vertex_partition 		= self.graph.nodes[self.fixed_vertex]["partition"]

			for neighbor in self.graph.neighbors(self.fixed_vertex):

				edge_weight 			= self.graph.edges[neighbor, self.fixed_vertex]["weight"]
				
				if fixed_vertex_partition != self.partitions[-1]:
					qubit 				= [qubit_map[neighbor][fixed_vertex_partition]]
					self.gate_i_z_1(gamma * edge_weight, qubit)

				else:
					for partition in self.partitions[:-1]:
						qubit 			= [qubit_map[neighbor][partition]]
						self.gate_i_z_1(- gamma * edge_weight, qubit)



			#-----------------------------------------------------------------------------------
			# Unitary operators associated to penalty terms of constraint violation (R-QUBO)
			#-----------------------------------------------------------------------------------
			for vertex in self.reduced_vertices:
				for (partition1, partition2) in itertools.combinations(self.partitions[:-1], 2):
					qubits 				= [qubit_map[vertex][partition1], qubit_map[vertex][partition2]]
					self.gate_i_z_2(gamma * self.penalty_coef[vertex], qubits)


		#---------------------------------------------------------------------------------------
		# Mixing operator (QUBO, PUBO, and R-QUBO)
		#---------------------------------------------------------------------------------------
		self.qaoa_circuit.rx(2 * beta, all_qubits)

	#--------------------------------------------------------------------------------------------
	# Measurement (QUBO, PUBO, and R-QUBO)
	#--------------------------------------------------------------------------------------------
	self.qaoa_circuit.measure(all_qubits, all_qubits)

	noise_model 		= NoiseModel()
	error_u1 			= depolarizing_error(self.Params.QAOA_Gates_Noise, 1)
	noise_model.add_all_qubit_quantum_error(error_u1, ['u1', 'rx', 'unitary'])

	backend 			= Aer.get_backend("qasm_simulator", fusion_max_qubit=3, precision="single") #qasm_simulator
	# backend 			= Aer.get_backend("qasm_simulator") #qasm_simulator

	shots				= self.Params.QAOA_Num_Shots

	
	simulate 			= execute(self.qaoa_circuit, backend=backend, shots=shots, noise_model=noise_model)
	qaoa_results		= simulate.result()
	
	self.counts 		= qaoa_results.get_counts()


	#--------------------------------------------------------------------------------------------
	# Calculate the average and the best objective values and the best solution 
	#--------------------------------------------------------------------------------------------
	self.cal_avg_best_sol(angles)
	self.qaoa_circuit_depth				= self.qaoa_circuit.depth()

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



