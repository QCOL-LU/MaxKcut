import numpy as np
from scipy.optimize import minimize, shgo, brute

import sys 
import time
from copy import deepcopy

from max_k_cut.quantum_methods._max_k_cut_qaoa_circuits import qaoa_expected_value
from max_k_cut.general_methods._plot_figures import plot_qaoa_solutions_dist
from max_k_cut.general_methods._print_methods import print_qaoa_results_summary

from max_k_cut.dfo_methods._brute import optimize_brute_search
from max_k_cut.dfo_methods._zo_svrg import optimize_zo_svrg

#================================================================================================
# Solve the max k-cut problem using QAOA
#================================================================================================
def solve_max_k_cut_qaoa(self):

	#--------------------------------------------------------------------------------------------
	# Parameters setting 
	#--------------------------------------------------------------------------------------------
	if not self.Params.QAOA_Angles:
		self.Params.QAOA_Angles		= [np.pi for _ in range(2*self.Params.QAOA_Num_Levels)]
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

	#--------------------------------------------------------------------------------------------
	# Optimize parameters of QAOA (gamma and beta angles)
	#--------------------------------------------------------------------------------------------
	if self.Params.QAOA_Optimize == True:

		self.qaoa_iter 				= 0

		if self.Params.QAOA_Scipy_Optimizer in ["L-BFGS-B", "TNC", "SLSQP", "Powell", "trust-constr"]:

			bounds 					= tuple([(0, 2 * np.pi) for _ in range(2*self.Params.QAOA_Num_Levels)])
			optimizer_results 		= minimize(fun=self.qaoa_expected_value, 
												x0=self.Params.QAOA_Angles, 
												method=self.Params.QAOA_Scipy_Optimizer, 
												bounds=bounds, 
												tol=self.Params.QAOA_Opt_Tol)
			locally_optimal_angles 	= optimizer_results.x

		elif self.Params.QAOA_Scipy_Optimizer == "Nelder-Mead":

			bounds 					= tuple([(0, 2 * np.pi) for _ in range(2*self.Params.QAOA_Num_Levels)])
			optimizer_results 		= minimize(fun=self.qaoa_expected_value, 
												x0=self.Params.QAOA_Angles, 
												method=self.Params.QAOA_Scipy_Optimizer, 
												tol=self.Params.QAOA_Opt_Tol)
			locally_optimal_angles 	= optimizer_results.x

		elif self.Params.QAOA_Scipy_Optimizer == "shgo":
			bounds 					= tuple([(0, 2 * np.pi) for _ in range(2*  self.Params.QAOA_Num_Levels)])
			optimizer_results 		= shgo(func=self.qaoa_expected_value, 
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

			optimizer_results 		= minimize(fun=self.qaoa_expected_value, 
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
		self.qaoa_expected_value(self.Params.QAOA_Angles)
		

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
