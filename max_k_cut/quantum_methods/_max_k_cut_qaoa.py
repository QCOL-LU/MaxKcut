import numpy as np
from scipy.optimize import minimize, shgo, brute

import sys 
import time
from copy import deepcopy

from max_k_cut.quantum_methods._max_k_cut_qaoa_dependencies import cal_obj_from_sol, convert_string_sol_to_sorted_sol, make_sol_feasible, cal_avg_best_sol, cal_avg_best_sol_feasible, gate_i_zz, gate_i_z_1, gate_i_z_2
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, execute, Aer, transpile, IBMQ

from max_k_cut.quantum_methods._max_k_cut_qaoa_circuits import qaoa_expected_value
from max_k_cut.general_methods._plot_figures import plot_qaoa_solutions_dist
from max_k_cut.general_methods._print_methods import print_qaoa_results_summary

from max_k_cut.dfo_methods._brute import optimize_brute_search
from max_k_cut.dfo_methods._zo_svrg import optimize_zo_svrg
from qiskit_ibm_runtime import QiskitRuntimeService, Estimator, Session, Sampler, Options

from qiskit.quantum_info import SparsePauliOp
import logging

#================================================================================================
# Solve the max k-cut problem using QAOA
#================================================================================================
def solve_max_k_cut_qaoa(self):

	#--------------------------------------------------------------------------------------------
	# Parameters setting 
	#--------------------------------------------------------------------------------------------
	if not self.Params.QAOA_Angles:
		self.Params.QAOA_Angles		= [0 for _ in range(2*self.Params.QAOA_Num_Levels)]
	else:
		self.Params.QAOA_Num_Levels = int(len(self.Params.QAOA_Angles)/2)

	self.reduced_vertices 	= [vertex for vertex in self.vertices if (vertex != self.fixed_vertex)]

	self.qaoa_best_obj_value		= - sys.maxsize
	self.qaoa_feasible_best_obj_value		= - sys.maxsize

	self.qaoa_skewness 				= 0.
	self.qaoa_best_skewness 		= 0.0

	self.qoao_opt_start_time		= time.time()
	self.qoao_opt_end_time 			= time.time()

	self.seed 						= 0

	if self.Params.Is_Simulator == True:
		backend 					= Aer.get_backend("qasm_simulator", fusion_max_qubit=3, precision="single")
		job 						= execute(self.qaoa_circuit, backend=backend, shots=shots, noise_model=noise_model)
	# else:
	# # 	# if self.Params.QLSA_First_Run == True:
	# 	IBMQ.save_account(self.Params.Token, overwrite=True)
	# 	IBMQ.load_account()
	# 	provider 					= IBMQ.get_provider(hub="ibm-q", group="open", project="main") 
	# 	backend 					= provider.get_backend(self.Params.QC_Name)

	shots				= self.Params.QAOA_Num_Shots


	self.qaoa_expected_value()
	parameter_values  				= {str(param): ind for ind, param in enumerate(self.qaoa_circuit.parameters)}


	gammas_keys 					= {level: 'gamma_{}'.format(level) for level in range(self.Params.QAOA_Num_Levels)}
	betas_keys 						= {level: 'beta_{}'.format(level) for level in range(self.Params.QAOA_Num_Levels)}

	self.qaoa_circuit_depth			= self.qaoa_circuit.depth()
	

	service 						= QiskitRuntimeService(channel="ibm_quantum", token=self.Params.Token)
	options 						= Options()
	options.execution.shots 		= self.Params.QAOA_Num_Shots

	self.error_iter 						= 0

	
	
	

	
	
	with Session(service=service, backend=self.Params.QC_Name) as session:
		sampler 					= Sampler(session=session)

		
		self.First_run						= True
		

		def evaluate_expectation(angles):
			
			self.gammas 						= np.array([angles[parameter_values[gammas_keys[level]]] for level in range(self.Params.QAOA_Num_Levels)])
			self.betas 							= np.array([angles[parameter_values[betas_keys[level]]] for level in range(self.Params.QAOA_Num_Levels)])

			temp_start_time 					= time.time()


			while self.error_iter < 10:
				try:
					job 						= sampler.run(circuits=self.qaoa_circuit, parameter_values=[list(angles)])
					qaoa_results 				= job.result()
				
				except:
					self.print_error_message()
					self.error_iter 					+= 1

				
				else:
					
					self.counts 				= qaoa_results.quasi_dists[0].binary_probabilities()
					self.cal_avg_best_sol(np.append(self.gammas, self.betas))
				
					if self.First_run: 
						self.qoao_opt_start_time		= time.time() 
						self.First_run 					= False
						self.waiting_time 				= self.qoao_opt_start_time - temp_start_time
						

					#------------------------------------------------------------------------------------
					# Print the results of QAOA optimizer at each iteration  
					#------------------------------------------------------------------------------------
					if self.Params.QAOA_Optimize == True:
						self.print_qaoa_optimizer_iter()	

					else:
						self.qoao_opt_end_time 			= time.time()
						self.qaoa_opt_total_time 		= self.qoao_opt_end_time - self.qoao_opt_start_time

					self.error_iter 							= 0
					break
						

			return -self.qaoa_avg_obj_value


		#--------------------------------------------------------------------------------------------
		# Optimize parameters of QAOA (gamma and beta angles)
		#--------------------------------------------------------------------------------------------
		if self.Params.QAOA_Optimize == True:

			self.qaoa_iter 				= 0
			bounds 						= tuple([(0, 2 * np.pi) for _ in range(2*  self.Params.QAOA_Num_Levels)])

			if self.Params.QAOA_Scipy_Optimizer in ["L-BFGS-B", "TNC", "SLSQP", "Powell", "trust-constr"]:

				optimizer_results 		= minimize(fun=evaluate_expectation, 
													x0=self.Params.QAOA_Angles, 
													method=self.Params.QAOA_Scipy_Optimizer, 
													bounds=bounds, 
													tol=self.Params.QAOA_Opt_Tol)
				locally_optimal_angles 	= optimizer_results.x

			elif self.Params.QAOA_Scipy_Optimizer == "Nelder-Mead":

				optimizer_results 		= minimize(fun=evaluate_expectation, 
													x0=self.Params.QAOA_Angles, 
													method=self.Params.QAOA_Scipy_Optimizer, 
													tol=self.Params.QAOA_Opt_Tol)
				locally_optimal_angles 	= optimizer_results.x

			elif self.Params.QAOA_Scipy_Optimizer == "shgo":
				
				optimizer_results 		= shgo(func=evaluate_expectation, 
												bounds=bounds,
												options={"f_tol": self.Params.QAOA_Opt_Tol})
				locally_optimal_angles 	= optimizer_results.x

			elif self.Params.QAOA_Scipy_Optimizer == "brute":
				optimizer_results 		= optimize_brute_search(self)

				self.plot_qaoa_level_one(optimizer_results.gammas, optimizer_results.betas, optimizer_results.objectives)
				locally_optimal_angles 	= optimizer_results.locally_optimal_angles


			elif self.Params.QAOA_Scipy_Optimizer == "zo_svrg":
				optimizer_results 		= optimize_zo_svrg(self)

			else:

				optimizer_results 		= minimize(fun=evaluate_expectation, 
													x0=self.Params.QAOA_Angles, 
													method=self.Params.QAOA_Scipy_Optimizer, 
													tol=self.Params.QAOA_Opt_Tol)
			
				locally_optimal_angles 	= optimizer_results.x

			with open(self.filename, 'a') as file:
				if self.Params.QAOA_Scipy_Optimizer != "brute":
					self.my_print(file, "\n" + optimizer_results.message + " It stopped after "+ str(self.qaoa_iter) + " iterations.\n")
				else:
					self.my_print(file, optimizer_results.message)
		else: 
			evaluate_expectation(self.Params.QAOA_Angles)
			

		#--------------------------------------------------------------------------------------------
		# Improve the resulting binary solution of QAOA and make it feasible in case of infeasibility
		#--------------------------------------------------------------------------------------------
		self.cal_avg_best_sol_feasible()
		

		#--------------------------------------------------------------------------------------------
		# Print the summary of results of QAOA 
		#--------------------------------------------------------------------------------------------
		if self.Params.QAOA_Verbosity > 0:

			self.plot_qaoa_solutions_dist(feasible_sol=False, num_bins=20)
			self.plot_qaoa_solutions_dist(feasible_sol=True, num_bins=10)

			self.print_qaoa_results_summary()

		#--------------------------------------------------------------------------------------------
		# Update the selected partitions for each of vertices in the original graph
		#--------------------------------------------------------------------------------------------

		best_solution  		= deepcopy(self.qaoa_feasible_best_solution) if self.qaoa_feasible_best_obj_value > self.qaoa_best_obj_value else deepcopy(self.qaoa_best_solution)
		
		for vertex in self.vertices:
			for partition in self.partitions:
				if best_solution[vertex][partition] == 1:
					self.graph.nodes[vertex]["partition"] 	= partition

