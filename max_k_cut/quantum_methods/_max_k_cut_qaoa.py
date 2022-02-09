import numpy as np
from scipy.optimize import minimize, shgo, brute

import sys 
import time
from copy import deepcopy

from quantum_methods._max_k_cut_qaoa_circuits import qaoa_expected_value
from general_methods._plot_figures import plot_qaoa_solutions_dist
from general_methods._print_methods import print_qaoa_results_summary

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
	self.qoao_opt_start_time		= time.time()
	self.qoao_opt_end_time 			= time.time()

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
			bounds 					= [(0,  np.pi) for _ in range(2 * self.Params.QAOA_Num_Levels)]
			# bounds					= [(.5, .9), (.6, 1), (1.8, 2.2), (4.8, 5.2)]
			Ns						= self.Params.QAOA_Brute_Num_Samples

			# optimizer_results 		= brute(func=self.qaoa_expected_value, 
			# 								ranges=bounds,
			# 								Ns=10)

			iter_finish				= 100

			if self.Params.QAOA_Num_Levels == 1:
				gammas, betas		= np.meshgrid(np.linspace(bounds[0][0]/np.pi, bounds[0][1]/np.pi, Ns), np.linspace(bounds[self.Params.QAOA_Num_Levels ][0]/np.pi, bounds[self.Params.QAOA_Num_Levels ][1]/np.pi, Ns))
				
				objectives 			= np.zeros((Ns, Ns))


			itr 					= 0
			total_itr 				= 0
			for g_ind, gamma1 in enumerate(np.linspace(*bounds[0], Ns) ):
				for b_ind, beta1 in enumerate(np.linspace(*bounds[self.Params.QAOA_Num_Levels ], Ns)):

					if self.Params.QAOA_Num_Levels == 1:
						old_best_exp_obj			= self.qaoa_best_avg_obj_value

						arg 						= (gamma1, beta1)
						obj 						= - self.qaoa_expected_value(arg)
						total_itr 					+= 1
						objectives[g_ind][b_ind]	= obj

						if self.qaoa_best_avg_obj_value != old_best_exp_obj:
							itr 					= 0
						else:
							itr 					+= 1

					if itr > iter_finish: break


					if self.Params.QAOA_Num_Levels >= 2:
						for gamma2 in np.linspace(*bounds[1], Ns):
							for beta2 in np.linspace(*bounds[self.Params.QAOA_Num_Levels + 1], Ns):

								if self.Params.QAOA_Num_Levels == 2:
									arg 	= (gamma1, gamma2, beta1, beta2)
									obj 	= - self.qaoa_expected_value(arg)

								if self.Params.QAOA_Num_Levels == 3:
									for gamma3 in np.linspace(*bounds[2], Ns):
										for beta3 in np.linspace(*bounds[self.Params.QAOA_Num_Levels + 2], Ns):

										
											arg 	= (gamma1, gamma2, gamma3, beta1, beta2, beta3)
											obj 	= - self.qaoa_expected_value(arg)

					
				if itr > iter_finish: break
				

			self.plot_qaoa_level_one(gammas, betas, objectives)
			locally_optimal_angles 		= self.best_angles

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
				self.my_print(file, "The brute search stopped after "+ str(total_itr) + " iterations.")
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
