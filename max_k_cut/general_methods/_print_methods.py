import networkx as nx 
import numpy as np

import time
import os

import matplotlib.pyplot as plt 
from copy import deepcopy

from max_k_cut.main.ParametersDefault import Parameters
import sys


#================================================================================================
# Print updated paramters
#================================================================================================
def my_print(self, file, text="", sep=' ', end='\n'):
	original_stdout 	= sys.stdout
	sys.stdout 			= file
	print(text)
	sys.stdout 			= original_stdout
	print(text)

def print_paramters(self):
	self.update_parameters()

	default_paramters 	= Parameters()
	with open(self.filename, 'a') as file:
		self.my_print(file, 50*"=")
		self.my_print(file, "Instance parameters")
		self.my_print(file, 50*"-")

		self.my_print(file, "{:<30}{:.0f}".format('Number of vertices:', self.num_vertices) )
		self.my_print(file, "{:<30}{:.0f}".format('Number of edges:', self.num_edges) )
		self.my_print(file, "{:<30}{:.0f}".format('Number of partitions:', self.num_partitions) )
		self.my_print(file)
		self.my_print(file, "{:<30}{:.0f}".format('Density (%):', self.density) )
		self.my_print(file, "{:<30}{:}".format('Name of instance:', self.name_specifier) )
		self.my_print(file)
		self.my_print(file, "{:<30}{:}".format('Is planar:', self.is_planar) )
		self.my_print(file, "{:<30}{:}".format('Is Chordal:', self.is_chordal ) )
		self.my_print(file)
		self.my_print(file, "{:<30}{:.0f}".format('Triangles Density (%):', self.triangles_density) )
		self.my_print(file, "{:<30}{:.0f}".format('Min_Maximal_Matching (%):', self.min_maximal_matching) )
		self.my_print(file, "{:<30}{:.2f}".format('Global Efficiency:', self.global_efficiency) )
		self.my_print(file)
		self.my_print(file, "{:<30}{:.0f}".format('Core Number:', self.core_number) )
		self.my_print(file, "{:<30}{:.0f}".format('Largest Component Size:', self.largest_component) )
		self.my_print(file)
		self.my_print(file, "{:<30}{:.0f}".format('Cut Vertices Number:', self.num_cut_vertices) )
		self.my_print(file, "{:<30}{:.0f}".format('Largest BiComponent Size:', self.largest_biconnected_component) )
		
		
		self.my_print(file, 50*"=")
		self.my_print(file)

		self.my_print(file, 50*"=")
		self.my_print(file, "Solver parameters")
		self.my_print(file, 50*"-")

		self.my_print(file, "{:30}{}".format("Method:", self.Params.Method) )

		self.my_print(file, "{:30}{:}".format("Peel:", self.Params.Peel) )

		self.my_print(file, "{:30}{:}".format("Decompose:", self.Params.Decompose) )

		self.my_print(file, "{:30}{:}".format("Fold:", self.Params.Fold) )			

		if (default_paramters.Relaxed != self.Params.Relaxed):
			self.my_print(file, "{:30}{:}".format("Relaxed:", self.Params.Relaxed) )

		

		if self.Params.Method == "BQO":

			self.my_print(file, "{:30}{:}".format("Rounding Heuristic:", self.Params.Rounding_Heuristic) )

			self.my_print(file, "{:30}{:}".format("Curvature Type:", self.Params.Curvature_Type) )

			self.my_print(file, "{:30}{:}".format("Curvature Method:", self.Params.Curvature_Method) )

			self.my_print(file, "{:30}{:}".format("Symmetry Breaking:", self.Params.Symmetry_Breaking) )


		if self.Params.Method == "A-MILO":

			self.my_print(file, "{:30}{:}".format("Clique Constraints:", self.Params.Clique_Constraints) )

			self.my_print(file, "{:30}{:}".format("Wheel Constraints:", self.Params.Wheel_Constraints) )

			self.my_print(file, "{:30}{:}".format("Bicycle Wheel Constraints:", self.Params.Bicycle_Wheel_Constraints) )


		elif self.Params.Method in ["QUBO", "R-QUBO", "PUBO", "C-QUBO", "CR-QUBO"]:
			self.my_print(file)

			self.my_print(file, "{:30}{:}".format("Quantum Computer:", self.Params.QC_Name) )

			self.my_print(file, "{:30}{:}".format("Naive Penalty Coef:", self.Params.Naive) )
			self.my_print(file, "{:30}{:}".format("Adjusted Penalty Coef:", self.Params.Adjusted_Penalty) )
			self.my_print(file, "{:30}{:}".format("Penalty Multiplier:", self.Params.Penalty_Increase) )
			

			self.my_print(file)
			
			self.my_print(file, "{:30}{:}".format("QAOA Num Levels:", self.Params.QAOA_Num_Levels) )
			self.my_print(file, "{:30}{:}".format("QAOA Num Shots:", self.Params.QAOA_Num_Shots) )
			self.my_print(file, "{:30}{:}".format("Gates Error Probability:", self.Params.QAOA_Gates_Noise) )

				
			
			self.my_print(file)
			self.my_print(file, "{:30}{:}".format("QAOA Optimize:", self.Params.QAOA_Optimize) )

			self.my_print(file, "{:30}{:}".format("QAOA Scipy Optimizer:", self.Params.QAOA_Scipy_Optimizer) )

			self.my_print(file, "{:30}{:}".format("QAOA Optimizer Tol:", self.Params.QAOA_Opt_Tol) )

			self.my_print(file, "{:30}{:}".format("QAOA Angles:", self.Params.QAOA_Angles) )

			

		if self.Params.Method not in ["QUBO", "R-QUBO", "PUBO", "MISDO", "MISDO2"]:

			self.my_print(file, "{:30}{:.0f}".format("Gurobi Time Limit:", self.Params.Gurobi_TimeLimit) )

			self.my_print(file, "{:30}{:.0f}".format("Gurobi MIPGap:", self.Params.Gurobi_MIPGap) )

			self.my_print(file, "{:30}{:.0f}".format("Gurobi Threads:", self.Params.Gurobi_Threads) )

			self.my_print(file, "{:30}{:.0f}".format("Gurobi Heuristics:", self.Params.Gurobi_Heuristics) )

			self.my_print(file, "{:30}{:.0f}".format("Gurobi Presolve:", self.Params.Gurobi_Presolve) )

			self.my_print(file, "{:30}{:.0f}".format("Gurobi Symmetry:", self.Params.Gurobi_Symmetry) )

			self.my_print(file, "{:30}{:.0f}".format("Gurobi Cuts:", self.Params.Gurobi_Cuts) )


		self.my_print(file, 50*"=")
		self.my_print(file)


#================================================================================================
# Print the summary of results of each instance
#================================================================================================
def print_results_summary(self):

	if self.Params.Verbosity > 0 or self.parent == None:
			
		self.objective_value 				= sum([self.graph.edges[edge]["weight"] for edge in self.edges if self.graph.nodes[edge[0]]["partition"] != self.graph.nodes[edge[1]]["partition"]]) 

		with open(self.filename, 'a') as file:
			self.my_print(file, 50*"=")
			self.my_print(file, "Summary of results of {:}".format(self.tree_node_name.replace("\n", "")))
			self.my_print(file, 50*"-")

			self.my_print(file, "{:<30}{:.0f}".format('Number of vertices:', self.num_vertices) )
			self.my_print(file, "{:<30}{:.0f}".format('Number of edges:', self.num_edges) )
			
			if self.Params.Method in ["RP-MILO", "C-RP-MILO"] and self.applied_operation == "solve":
				self.my_print(file, "{:<30}{:.0f}".format('Density of ext Chordal (%):', self.density_chordal_graph) )

			if self.applied_operation == "decompose":
				self.my_print(file, "{:<30}{:.0f}".format('Num of bi-comp:', self.num_biconnected_component) )

			self.my_print(file)

			self.my_print(file, "{:<30}{:.0f}".format('Vertex num in largest comp:', self.largest_subgraph[0]) )
			self.my_print(file, "{:<30}{:.0f}".format('Edge num in largest comp:', self.largest_subgraph[1]) )

			self.my_print(file)

			self.my_print(file, "{:<30}{:.2f}".format('Pre-processing running time:', self.preprocessing_total_time) )
			self.my_print(file, "{:<30}{:.2f}".format('Total solver time:', self.total_solver_time) )
			self.my_print(file, "{:<30}{:.2f}".format('Running time:', self.end_total_run_time - self.start_total_run_time) )

			self.my_print(file)
			
			self.my_print(file, "{:<30}{:.2f}".format('Upper bound: ', self.upper_bound) )
			self.my_print(file, "{:<30}{:.2f}".format('Objective value: ', self.objective_value) )
			

			if self.Params.Verbosity == 2:
				self.my_print(file)
				for partition in self.partitions:
					assigned_vertices 			= [vertex for vertex in self.vertices if self.graph.nodes[vertex]["partition"] == partition]
					assigned_vertices_string 	= (str(assigned_vertices)[1:-1] if bool(assigned_vertices) == True else '')
					self.my_print(file, "Partition {:}: {:<16} {:<20}".format(partition," " ,assigned_vertices_string) )

			self.my_print(file, 50*"=")
			if self.parent != None: self.my_print(file)





#================================================================================================
# Print the summary of results of QAOA 
#================================================================================================
def print_qaoa_results_summary(self):
	if self.Params.Verbosity > 0 or self.parent == None:
		modified_qaoa_best_sol							= self.make_sol_feasible(self.qaoa_best_solution)
		temp, modified_qaoa_best_sol_obj_value, penalty, const_viol 		= self.cal_obj_from_sol(modified_qaoa_best_sol)
			
		with open(self.filename, 'a') as file:
			self.my_print(file)
			self.my_print(file, 50*"=")
			self.my_print(file, "Summary of results QAOA")
			self.my_print(file, 50*"-")

			bqo_obj 			= self.qaoa_best_avg_total_penalty + self.qaoa_best_avg_obj_value

			self.my_print(file, "{:<30}{:.2f}".format("Avg of QAOA obj (BQO obj):", bqo_obj))
			self.my_print(file, "{:<30}{:.2f}".format("Avg of QAOA obj (penalty):", - self.qaoa_best_avg_total_penalty))
			self.my_print(file, "{:<30}{:.2f}".format("Avg of constraint violation:", self.qaoa_best_avg_constraint_violation))
			self.my_print(file)
 
			self.my_print(file, "{:<30}{:.2f}".format("Feasible percentage (%):", self.qaoa_best_pure_feasible_percentage))
			if self.qaoa_best_avg_pure_feasible_obj_value == "NA":
				self.my_print(file, "{:<30}{:}".format("Avg of pure feasible obj:", self.qaoa_best_avg_pure_feasible_obj_value))

			else:
				self.my_print(file, "{:<30}{:.2f}".format("Avg of pure feasible obj:", self.qaoa_best_avg_pure_feasible_obj_value))
			self.my_print(file)
			self.my_print(file, "{:<30}{:.2f}".format("Avg of tight R-QUBO QAOA obj (penalty):", self.qaoa_best_avg_rqubo_total_penalty))
			self.my_print(file, "{:<30}{:.2f}".format("Avg of tight R-QUBO QAOA obj:", self.qaoa_best_avg_rqubo_tight_obj_value))
			self.my_print(file)
			self.my_print(file, "{:<30}{:.2f}".format("Avg of tight QAOA obj (penalty):",self.qaoa_best_avg_tight_obj_value - bqo_obj))
			self.my_print(file, "{:<30}{:.2f}".format("Avg of tight QAOA obj:", self.qaoa_best_avg_tight_obj_value))
			self.my_print(file)
			self.my_print(file, "{:<30}{:.2f}".format("Avg of QAOA obj:", self.qaoa_best_avg_obj_value))
			
			self.my_print(file, "{:<30}{:.2f}".format("STD of QAOA obj:", self.qaoa_std_obj_value))

			self.my_print(file)
			self.my_print(file, "{:<30}{:.2f}".format("Avg of QAOA obj (feasible):", self.qaoa_feasible_best_avg_obj_value))
			self.my_print(file, "{:<30}{:.2f}".format("STD of QAOA obj (feasible):", self.qaoa_feasible_std_obj_value))

			if self.Params.Verbosity == 2:
				self.my_print(file)
				self.my_print(file, "{:<30}{:}".format("Gamma angles:", ', '.join(('%.2f'% f) for f in self.best_angles[:self.Params.QAOA_Num_Levels]) ) )
				self.my_print(file, "{:<30}{:}".format("Beta angles:", ', '.join(('%.2f'% f) for f in self.best_angles[self.Params.QAOA_Num_Levels:]) ) )
				
				
			self.my_print(file)
			self.my_print(file, "{:<30}{:.2f}".format("Best QAOA obj:", self.qaoa_best_obj_value))
			self.my_print(file, "{:<30}{:.2f}".format("Modified best QAOA obj:", modified_qaoa_best_sol_obj_value))

			self.my_print(file, "{:<30}{:.2f}".format("Best QAOA obj (feasible):", self.qaoa_feasible_best_obj_value))
			self.my_print(file)
			
			self.my_print(file, "{:<30}{:.2f}s".format("QAOA wating time (s):", self.waiting_time))
			self.my_print(file, "{:<30}{:.2f}".format("QAOA total time (s):", self.qaoa_opt_total_time))
			self.my_print(file, "{:<30}{:.0f}".format("QAOA circuit depth:", self.qaoa_circuit_depth))

			self.my_print(file, 50*"=")
			self.my_print(file)




#================================================================================================
# Print the results of QAOA optimizer at each iteration 
#================================================================================================
def print_qaoa_optimizer_iter(self):
	np.set_printoptions(precision=2, formatter={'float': ('{:.'+str(2)+'f} ').format})

	if self.Params.Verbosity > 0 or self.parent == None:
		if self.qaoa_iter == 0:
			self.qaoa_temp_best 		= - sys.maxsize
			with open(self.filename, 'a') as file:
				self.my_print(file, " {:>6}  {:^20}{:^15}  {:^15}  {:^15}  {:^15}  {:^15}  {:^15} {:>5}     {:^15}  {:^15}".format( "iter", "avg-obj", "best avg-obj", "best-obj", "tight-obj", "rqubo-obj", "std", "skewness" ,"time", "gamma", "beta"))
		
		do_print_iter		= 	(self.Params.QAOA_Verbosity == 2) or \
								(self.Params.QAOA_Verbosity == 1 and \
									time.time() - self.qoao_opt_end_time >= self.Params.QAOA_Opt_Print_Time)

		self.qaoa_opt_total_time 		= time.time() - self.qoao_opt_start_time

		
		if do_print_iter == True:
			self.qoao_opt_end_time 		= time.time()

			if self.qaoa_temp_best < self.qaoa_avg_obj_value: 
				self.qaoa_temp_best 	= self.qaoa_avg_obj_value

				star 		= "*" 
			else:
				star 		= " "

			
			with open(self.filename, 'a') as file:
				self.my_print(file, "{:}{:>6d}  {:>15.5e}  {:>15.5e}  {:>15.5e}  {:>15.5e}  {:>15.5e}  {:>15.5e}  {:>15.5e}  {:>5.0f}s      {:^15}  {:^15}".format(star, self.qaoa_iter, self.qaoa_avg_obj_value, self.qaoa_best_avg_obj_value,  self.qaoa_best_obj_value, self.qaoa_avg_tight_obj_value, self.qaoa_avg_rqubo_tight_obj_value, self.qaoa_std_obj_value, self.qaoa_skewness,self.qaoa_opt_total_time, str(self.gammas)[1:-1], str(self.betas)[1:-1]))

	self.qaoa_iter 					= self.qaoa_iter + 1

#================================================================================================
# Print the results of QAOA optimizer at each iteration 
#================================================================================================
def print_gurobi_results_summary(self):
	if self.Params.Verbosity > 0 or self.parent == None:
		with open(self.filename, 'a') as file:
			self.my_print(file, 50*"=")
			self.my_print(file, "Summary of results of {:} model of {:}".format(self.Params.Method, self.tree_node_name.replace("\n", "")))
			self.my_print(file, 50*"-")

			self.my_print(file, "{:<30}{:}".format('Model status: ', self.gurobi_model_status) )
			self.my_print(file, "{:<30}{:.0f}".format('Explored B&B nodes:',  self.gurobi_BB_nodes) )
			self.my_print(file, "{:<30}{:.1f}".format('Solver running time:', self.gurobi_end_time - self.gurobi_start_time) )
			
			self.my_print(file)
			self.my_print(file, "{:<30}{:.2f}".format('Optimality gap (%): ', 100*self.gurobi_MIPGap) )	
			self.my_print(file, "{:<30}{:.2f}".format('Upper bound: ', self.gurobi_ObjBound) )	
			self.my_print(file, "{:<30}{:.2f}".format('Objective value: ', self.gurobi_obj_value) )	

			if self.Params.Rounding_Heuristic == True and self.Params.Method == "BQO" and self.Params.Relaxed == True:
				heuristic_obj 		= sum([1 for edge in self.graph.edges if self.graph.nodes[edge[0]]["partition"] != self.graph.nodes[edge[1]]["partition"]])
				self.my_print(file, "{:<30}{:.0f}".format('Heuristic obj value: ', heuristic_obj) )

			

			self.my_print(file, 50*"=")
			self.my_print(file)



def print_error_message(self):
	with open(self.filename, 'a') as file:
		self.my_print(file)
		self.my_print(file, "An error occurred in the IBM quantum computer!")
		self.my_print(file)




