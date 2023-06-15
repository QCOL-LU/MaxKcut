from max_k_cut.main.ParametersDefault import Parameters
import networkx as nx 

import matplotlib.pyplot as plt 
import time
import gurobipy as gp
from math import comb

import os
import sys



from ..formulations import *
from ..quantum_methods import *
from ..decomposition_methods import *
from ..general_methods import *



from max_k_cut.dfo_methods._cal_gradient import *




#====================================================================================================
# The class of instance; it contains different methods
#====================================================================================================
class Instance(Parameters):
	
	#================================================================================================
	# Class initialization
	#================================================================================================
	def __init__(self, graph, fixed_vertex_partition=(None, None), parameters=None, name_specifier="", parent=None, indentifier=0):

		self.parent 				= parent
		self.graph 					= graph 				# networkx weighted graph
		
		self.indentifier 			= indentifier
		self.decomposition_tree		= nx.DiGraph()

		self.edges					= [tuple(sorted(edge)) for edge in list(graph.edges)]
		self.num_edges				= len(self.edges)

		self.vertices 				= sorted(list(graph.nodes) )
		self.num_vertices 			= len(self.vertices)


		self.objective_value 		= sys.maxsize
		self.upper_bound			= 0
		self.total_solver_time		= 0
		
		self.preprocessing_total_time 				= 0
		self.qaoa_best_avg_obj_value 				= - sys.maxsize
		self.qaoa_feasible_best_avg_obj_value 		= - sys.maxsize


		
		self.largest_subgraph		= (0, 0)								# (num_vertices, num_edges)

		if nx.is_weighted(self.graph) == False: 
			for edge in self.edges:
				self.graph.edges[edge]["weight"] 	= 1

	
		self.total_weights 			= sum([self.graph.edges[edge]["weight"] for edge in self.edges]) 

		for vertex in self.vertices:
			neighbor_weights 						= [self.graph.edges[incident_edge]["weight"] for incident_edge in self.graph.edges(vertex)]
			self.graph.nodes[vertex]["pos-weight"]	= sum([weight for weight in neighbor_weights if weight > 0])
			self.graph.nodes[vertex]["neg-weight"]	= sum([weight for weight in neighbor_weights if weight < 0])

		
		nx.set_node_attributes(self.graph, -1, "partition")



		self.num_traversed_graphs 		= self.indentifier 
		self.applied_operation 			= None  	# 1: solve    2: peel    3: decompose    4: fold 	

		if self.num_vertices > 0:
			self.fixed_vertex  			= self.vertices[0] if fixed_vertex_partition[0] == None else fixed_vertex_partition[0]
		

		if self.parent == None:
			self.Params 				= Parameters() if parameters == None else parameters
			self.name_specifier  		= name_specifier

			self.density				= 200 * self.num_edges / ( (self.num_vertices - 1) * self.num_vertices ) if self.num_vertices > 1 else 0
			self.triangles_density		= int(sum(nx.triangles(self.graph).values()) / 3 / comb(self.num_vertices,3) * 100) 
			self.min_maximal_matching	= int(len(nx.algorithms.approximation.min_maximal_matching(self.graph)) / self.num_vertices * 200)
			self.global_efficiency 		= nx.global_efficiency(self.graph)
			
			self.is_planar 				= nx.check_planarity(self.graph)[0]
			self.is_chordal 			= nx.is_chordal(self.graph)

			self.core_number			= max(nx.core_number(self.graph).values() ) + 1 
			self.largest_component 		= len(max(nx.connected_components(self.graph), key=len))

			self.num_cut_vertices 				= sum(1 for _ in nx.articulation_points(self.graph)) 
			self.largest_biconnected_component	= len(max(nx.biconnected_components(self.graph), key=len))
			
		else:
			self.Params 				= self.parent.Params
			self.name 					= self.parent.name
			self.directory 				= self.parent.directory
			self.figure_path 			= self.parent.figure_path
			self.result_path 			= self.parent.result_path
			self.filename 				= self.parent.filename
			self.num_partitions 		= self.parent.num_partitions
			self.partitions 			= self.parent.partitions
			

			if self.num_vertices > 0:
				self.graph.nodes[self.fixed_vertex]["partition"]		= 0 if fixed_vertex_partition[1] == None else self.parent.graph.nodes[self.fixed_vertex]["partition"]


	#================================================================================================
	# Update parameters
	#================================================================================================	
	def update_parameters(self):
		self.num_partitions 	= self.Params.Num_Partitions
		self.partitions 		= [partition for partition in range(self.num_partitions)]

		self.num_partitions 	= self.Params.Num_Partitions
		self.partitions 		= [partition for partition in range(self.num_partitions)]

		self.category				= "d1" if self.density <= 5 else ("d2" if self.density <= 10 else ("d3" if self.density <= 25 else "d4")) + "_k" + str(self.num_partitions) 
		

		sizes    		= [49, 50, 100, 150, 200, 250, 251]
		densities 		= [5, 15, 30, 50, 70, 100]

		def find_position(my_list, number):
			if number <= my_list[0]: return my_list[0]

			for lb, ub in zip(my_list[:-1], my_list[1:]):
				if lb <= number < ub: return lb
			return my_list[-1]

		self.name  				= "k" + str(self.num_partitions).zfill(2)  
		self.directory 			= self.name + "_n" + str(find_position(sizes, self.num_vertices)).zfill(3) + "_d" + str(find_position(densities, int(self.density))).zfill(3)
		self.name 				= "k" + str(self.num_partitions).zfill(2) \
									+ "_n" + str(self.num_vertices).zfill(3) \
									+ "_d" + str(int(self.density)).zfill(3) + "_" + str(self.name_specifier)

		self.figure_path 		= "../figures/" + self.directory 
		self.result_path 		= "../results/" + self.directory

		self.filename 			= self.result_path + "/"+ self.name +"_" + self.Params.Method + ".txt"

		if not os.path.isdir(self.figure_path): os.mkdir(self.figure_path)
		if not os.path.isdir(self.result_path): os.mkdir(self.result_path)

		self.graph.nodes[self.fixed_vertex]["partition"]	= self.partitions[0] if self.parent == None else self.parent.graph.nodes[self.fixed_vertex]["partition"]


	#================================================================================================
	# Create the name of the associated node in the decomposition tree
	#================================================================================================
	def create_tree_node(self):
		mid_char 				= "-N" if self.applied_operation == None else "-" + self.applied_operation[0].capitalize()
		self.tree_node_name 	= str(self.indentifier) + mid_char + "\n" + str(self.num_vertices)
		self.decomposition_tree.add_node(self.tree_node_name) 


	#================================================================================================
	# Import functions 
	#================================================================================================

	#------------------------------------------------------------------------------------------------
	# Decomposition methods
	#------------------------------------------------------------------------------------------------
	from ..decomposition_methods._peel import peel, update_parent_peel
	from ..decomposition_methods._decompose import decompose, update_parent_decompose
	from ..decomposition_methods._fold import fold, update_parent_fold, folded_subgraph_solver
	
	from ..decomposition_methods._Lange_et_al_fold import lange_fold, update_parent_lange_fold

	from ..decomposition_methods._edge_based_fix import edge_based_fix, update_parent_edge_based_fix
	from ..decomposition_methods._Rehfeldt_et_al_fold import rehfeldt_fix, update_parent_rehfeldt_fix
	from ..decomposition_methods._twin_fix import twin_fix, update_parent_twin_fix

	#------------------------------------------------------------------------------------------------
	# Classical formulations
	#------------------------------------------------------------------------------------------------
	from ..formulations._max_k_cut_amilo import solve_max_k_cut_amilo
	from ..formulations._max_k_cut_pmilo import solve_max_k_cut_pmilo
	from ..formulations._max_k_cut_rpmilo import solve_max_k_cut_rpmilo
	from ..formulations._max_k_cut_crpmilo import solve_max_k_cut_crpmilo

	from ..formulations._max_k_cut_bqo import solve_max_k_cut_bqo
	from ..formulations._max_k_cut_rbqo import solve_max_k_cut_rbqo

	from ..formulations._max_k_cut_misdo import solve_max_k_cut_misdo
	from ..formulations._max_k_cut_misdo2 import solve_max_k_cut_misdo2

	from ..formulations._max_k_cut_qubo import solve_max_k_cut_qubo
	from ..formulations._max_k_cut_rqubo import solve_max_k_cut_rqubo

	from ..formulations._curvature_coefs import calculate_curvature_coefs

	#------------------------------------------------------------------------------------------------
	# QAOA method
	#------------------------------------------------------------------------------------------------
	from ..quantum_methods._max_k_cut_qaoa import solve_max_k_cut_qaoa
	from ..quantum_methods._max_k_cut_qaoa_circuits import qaoa_expected_value

	from ..quantum_methods._max_k_cut_qaoa_dependencies import cal_obj_from_sol, convert_string_sol_to_sorted_sol, make_sol_feasible, cal_avg_best_sol, cal_avg_best_sol_feasible, gate_i_zz, gate_i_z_1, gate_i_z_2
	from ..dfo_methods._cal_gradient import cal_gradient, cal_gradient_batch

	#------------------------------------------------------------------------------------------------
	# Plot methods
	#------------------------------------------------------------------------------------------------
	from ..general_methods._plot_figures import plot_qaoa_solutions_dist, plot_graph_problem, plot_qaoa_level_one
	from ..general_methods._plot_figures import plot_graph_solution, plot_decomposition_tree

	#------------------------------------------------------------------------------------------------
	# Print methods
	#------------------------------------------------------------------------------------------------
	from ..general_methods._print_methods import my_print, print_paramters, print_results_summary
	from ..general_methods._print_methods import print_qaoa_results_summary, print_qaoa_optimizer_iter
	from ..general_methods._print_methods import print_gurobi_results_summary


	
	#================================================================================================
	# Decompose and solve the problem
	#================================================================================================
	def solve(self):

		#--------------------------------------------------------------------------------------------
		# Method initialization
		#--------------------------------------------------------------------------------------------
		self.start_total_run_time 		= time.time()
		self.Params.Relaxed 			= True if self.Params.Rounding_Heuristic == True else self.Params.Relaxed

		if self.parent == None:			self.print_paramters()
		
		# is_peeling_allowed 				= self.Params.Peel and (True if self.parent == None else (not self.parent.applied_operation == "peel")) and (self.num_vertices > 0)
		# is_decompose_allowed  			= self.Params.Decompose and (True if self.parent == None else (not self.parent.applied_operation == "decompose")) and (self.num_vertices > 0)
		# is_folding_allowed  			= self.Params.Fold and (True if self.parent == None else (not self.parent.applied_operation == "fold")) and (self.num_vertices > 0)
		
		# is_fixing_twin_allowed  		= self.Params.Twin_Fix and (True if self.parent == None else (not self.parent.applied_operation == "twin-fix")) and (self.num_vertices > 0)
		# is_edge_based_fixing_allowed 	= self.Params.Edge_Based_Fix and (True if self.parent == None else (not self.parent.applied_operation == "edge-based-fix")) and (self.num_vertices > 0)
		# is_Rehfeldt_fixing_allowed 		= self.Params.Rehfeldt_Fix and (True if self.parent == None else (not self.parent.applied_operation == "Rehfeldt_fix")) and (self.num_vertices > 0)
		
		# is_Lange_folding_allowed 		= self.Params.Lange_Fold and (True if self.parent == None else (not self.parent.applied_operation == "Lange_fold")) and (self.num_vertices > 0)

		
		is_peeling_allowed 				= self.Params.Peel and (True if self.parent == None else  (self.num_vertices > 0) )
		is_decompose_allowed  			= self.Params.Decompose and (True if self.parent == None else (self.num_vertices > 0) )
		is_folding_allowed  			= self.Params.Fold and (True if self.parent == None else  (self.num_vertices > 0) )
		
		is_fixing_twin_allowed  		= self.Params.Twin_Fix and (True if self.parent == None else (self.num_vertices > 0) )
		is_edge_based_fixing_allowed 	= self.Params.Edge_Based_Fix and (True if self.parent == None else (self.num_vertices > 0) )
		is_Rehfeldt_fixing_allowed 		= self.Params.Rehfeldt_Fix and (True if self.parent == None else (self.num_vertices > 0) )
		
		is_Lange_folding_allowed 		= self.Params.Lange_Fold and (True if self.parent == None else  (self.num_vertices > 0) )


		is_graph_simplified 			= False 

		

		

		#--------------------------------------------------------------------------------------------
		# Apply peeling reducition algorithm if it is allowed
		#--------------------------------------------------------------------------------------------
		if is_peeling_allowed == True:

			is_graph_simplified			= self.peel()

			if is_graph_simplified == True:

				new_instance 			= Instance( graph=nx.Graph(self.graph.subgraph(self.k_core) ), 
													parent=self,
													indentifier=self.num_traversed_graphs + 1)
				self.create_tree_node()
				new_instance.solve()

		
		#--------------------------------------------------------------------------------------------
		# Apply biconnected decomposition algorithm if it is allowed
		#--------------------------------------------------------------------------------------------
		if not is_graph_simplified and is_decompose_allowed == True:

			is_graph_simplified 		= self.decompose()

			if is_graph_simplified == True:
				self.create_tree_node()

				for (ind, component) in enumerate(self.ordered_components):

					fixed_vertex		= self.fixed_vertices[ind]
					fixed_partition 	= self.graph.nodes[fixed_vertex]["partition"]

					new_instance 		= Instance( graph=nx.Graph(self.graph.subgraph(component)), 
													fixed_vertex_partition=(fixed_vertex, fixed_partition), 
													parent=self,
													indentifier=self.num_traversed_graphs + 1)

					new_instance.solve()

		#--------------------------------------------------------------------------------------------
		# Apply edge-based fixing algorithm if it is allowed
		#--------------------------------------------------------------------------------------------
		if not is_graph_simplified and is_edge_based_fixing_allowed == True:

			(is_graph_simplified, edge_based_fixed_graph) 	= self.edge_based_fix()


			if is_graph_simplified == True:
				new_instance 			= Instance( graph=edge_based_fixed_graph, 
													fixed_vertex_partition=(self.fixed_vertex, self.graph.nodes[self.fixed_vertex]["partition"]),
													parent=self, 
													indentifier=self.num_traversed_graphs + 1)
				self.create_tree_node()
				new_instance.solve()


		#--------------------------------------------------------------------------------------------
		# Apply Lange folding operation algorithm if it is allowed
		#--------------------------------------------------------------------------------------------
		if not is_graph_simplified and is_Lange_folding_allowed == True:

			(is_graph_simplified, Lange_folded_graph) = self.lange_fold()

			if is_graph_simplified == True:
				new_instance 			= Instance( graph=Lange_folded_graph, 
													fixed_vertex_partition=(self.fixed_vertex, self.graph.nodes[self.fixed_vertex]["partition"]),
													parent=self, 
													indentifier=self.num_traversed_graphs + 1)
				self.create_tree_node()
				new_instance.solve()

		#--------------------------------------------------------------------------------------------
		# Apply folding operation algorithm if it is allowed
		#--------------------------------------------------------------------------------------------
		if not is_graph_simplified and is_folding_allowed == True:

			(is_graph_simplified, folded_graph) 		= self.fold()

			if is_graph_simplified == True:
				new_instance 			= Instance( graph=folded_graph, 
													fixed_vertex_partition=(self.fixed_vertex, self.graph.nodes[self.fixed_vertex]["partition"]),
													parent=self, 
													indentifier=self.num_traversed_graphs + 1)
				self.create_tree_node()
				new_instance.solve()


				
		#--------------------------------------------------------------------------------------------
		# Apply Rehfeldt fixing operation algorithm if it is allowed
		#--------------------------------------------------------------------------------------------
		if not is_graph_simplified and is_Rehfeldt_fixing_allowed == True:

			(is_graph_simplified, Rehfeldt_fixed_graph) = self.rehfeldt_fix()

			if is_graph_simplified == True:
				new_instance 			= Instance( graph=Rehfeldt_fixed_graph, 
													fixed_vertex_partition=(self.fixed_vertex, self.graph.nodes[self.fixed_vertex]["partition"]),
													parent=self, 
													indentifier=self.num_traversed_graphs + 1)
				self.create_tree_node()
				new_instance.solve()


		


		#--------------------------------------------------------------------------------------------
		# Apply twin fixing algorithm if it is allowed
		#--------------------------------------------------------------------------------------------
		if not is_graph_simplified and is_fixing_twin_allowed == True:

			(is_graph_simplified, twin_fixed_graph) 	= self.twin_fix()


			if is_graph_simplified == True:
				new_instance 			= Instance( graph=twin_fixed_graph, 
													fixed_vertex_partition=(self.fixed_vertex, self.graph.nodes[self.fixed_vertex]["partition"]),
													parent=self, 
													indentifier=self.num_traversed_graphs + 1)
				self.create_tree_node()
				new_instance.solve()


		
		#--------------------------------------------------------------------------------------------
		# Solve the problem if the graph is not simplified
		#--------------------------------------------------------------------------------------------
		if is_graph_simplified == False and self.num_vertices > 0:


			if self.Params.Method == "BQO":								solver 		= solve_max_k_cut_bqo

			elif self.Params.Method == "R-BQO": 						solver 		= solve_max_k_cut_rbqo

			elif self.Params.Method == "A-MILO":						solver 		= solve_max_k_cut_amilo

			elif self.Params.Method == "P-MILO": 						solver 		= solve_max_k_cut_pmilo

			elif self.Params.Method == "RP-MILO": 						solver 		= solve_max_k_cut_rpmilo

			elif self.Params.Method == "C-RP-MILO": 					solver 		= solve_max_k_cut_crpmilo

			elif self.Params.Method == "MISDO":							solver 		= solve_max_k_cut_misdo

			elif self.Params.Method == "MISDO2":						solver 		= solve_max_k_cut_misdo2

			elif self.Params.Method == "C-QUBO": 						solver 		= solve_max_k_cut_qubo

			elif self.Params.Method == "CR-QUBO": 						solver 		= solve_max_k_cut_rqubo

			elif self.Params.Method in ["QUBO", "PUBO", "R-QUBO"]: 		solver 		= solve_max_k_cut_qaoa

				

			self.applied_operation 		= "solve"

			
			self.create_tree_node()

			solver(self)

			self.total_solver_time	= self.gurobi_end_time - self.gurobi_start_time if self.Params.Method not in ["QUBO", "PUBO", "R-QUBO"] else self.qoao_opt_end_time - self.qoao_opt_start_time
			

			self.upper_bound		= self.gurobi_ObjBound if self.Params.Method not in ["QUBO", "PUBO", "R-QUBO"] else sum(self.graph.edges[edge]["weight"] for edge in self.edges if self.graph.edges[edge]["weight"] > 0)
			self.largest_subgraph	= (self.num_vertices, self.num_edges)

		elif is_graph_simplified == False and self.num_vertices == 0:
			self.create_tree_node()

		
		#--------------------------------------------------------------------------------------------
		# Update the parent partitions (intermediate graphs)
		#--------------------------------------------------------------------------------------------
		if self.parent != None:
			
			if self.parent.applied_operation == "peel":				self.update_parent_peel()

			elif self.parent.applied_operation == "decompose":		self.update_parent_decompose()

			elif self.parent.applied_operation == "fold": 			self.update_parent_fold()

			elif self.parent.applied_operation == "lange-fold": 	self.update_parent_lange_fold()

			elif self.parent.applied_operation == "edge-based-fix": self.update_parent_edge_based_fix()

			elif self.parent.applied_operation == "rehfeldt-fix": 	self.update_parent_rehfeldt_fix()

			elif self.parent.applied_operation == "twin-fix": 		self.update_parent_twin_fix()

			

			self.parent.num_traversed_graphs 		= self.num_traversed_graphs

			self.decomposition_tree.add_edge(self.parent.tree_node_name, self.tree_node_name) 

			self.parent.decomposition_tree 			= nx.compose(self.parent.decomposition_tree, self.decomposition_tree)
		
			max_num_vertices						= max(self.parent.largest_subgraph[0], self.largest_subgraph[0])
			max_num_edges							= max(self.largest_subgraph[1], self.parent.largest_subgraph[1]) if self.largest_subgraph[0] == self.parent.largest_subgraph[0] else \
															(self.largest_subgraph[1] if max_num_vertices == self.largest_subgraph[0] else self.parent.largest_subgraph[1])

			self.parent.largest_subgraph			= (max_num_vertices, max_num_edges)

			self.parent.preprocessing_total_time 	+= self.preprocessing_total_time
			self.parent.total_solver_time			+= self.total_solver_time

		#--------------------------------------------------------------------------------------------
		# Plot the figures of the graph corresponding to the obtained solution and decomposition tree
		#--------------------------------------------------------------------------------------------
		else:
			
			self.plot_graph_solution()

			self.plot_decomposition_tree()
			

		self.end_total_run_time 				= time.time()
		self.print_results_summary()

	
		