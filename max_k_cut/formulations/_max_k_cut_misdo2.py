import networkx as nx 
import time
import threading
import itertools
import numpy as np
from copy import deepcopy

import sys
import mosek.fusion as mosek


#================================================================================================
# The MISDO formulation
#================================================================================================
def solve_max_k_cut_misdo2(self):

	self.gurobi_start_time 			= time.time()

	#--------------------------------------------------------------------------------------------
	# Model initialization
	#--------------------------------------------------------------------------------------------
	index_list 			= [ind for ind in range( self.num_vertices)]

	
	model 				= mosek.Model("MISDO2")
	Z 					= model.variable( mosek.Domain.inPSDCone( self.num_vertices) )

	neg_adjacency_matrix 	= np.zeros((self.num_vertices, self.num_vertices))
	objective_coefficient = (self.num_partitions - 1) / self.num_partitions

	for vertex1, vertex2 in self.edges:
		ind1, ind2 		= self.vertices.index(vertex1), self.vertices.index(vertex2)
		neg_adjacency_matrix[ind1, ind2] = - self.graph.edges[(vertex1, vertex2)]["weight"] * objective_coefficient

	objective 			= mosek.Expr.add(self.total_weights * objective_coefficient, mosek.Expr.sum(mosek.Expr.mulElm(Z, neg_adjacency_matrix) ) )


	model.objective(mosek.ObjectiveSense.Maximize, objective)


	for ind1, ind2 in itertools.combinations(index_list, 2):
		model.constraint(Z.index(ind1, ind2), mosek.Domain.greaterThan(-1 / (self.num_partitions - 1) ))

	for ind in index_list:
		model.constraint(Z.index(ind, ind), mosek.Domain.equalsTo(1))

	#----------------------------------------------------------------------------------------
	# Solve the SDO model
	#----------------------------------------------------------------------------------------
	time_limit 			= 60
	

	temp 				= sys.stdout
	sys.stdout 			= open(self.filename[:-4] + "_log.txt", "w")
	model.setLogHandler(sys.stdout)

	# 
	model.setSolverParam("numThreads", self.Params.Gurobi_Threads)
	model.setSolverParam("logFile", 1)
	model.setSolverParam("optimizerMaxTime", self.Params.Gurobi_TimeLimit)
	

	model.solve()
	sys.stdout.close() 
	sys.stdout 			= temp 
	
	#----------------------------------------------------------------------------------------
	# Extract the obtained solution
	#----------------------------------------------------------------------------------------
	if model.getPrimalSolutionStatus() == mosek.SolutionStatus.Feasible or model.getPrimalSolutionStatus() == mosek.SolutionStatus.Optimal:
		
		self.gurobi_model_status 	= 1
		self.gurobi_BB_nodes 		= 0 
		self.gurobi_MIPGap  		= 0 
		self.gurobi_ObjBound 		= model.dualObjValue()
		self.gurobi_obj_value 		= model.primalObjValue()

		z_star_edges 	= [(self.vertices[ind1], self.vertices[ind2]) for (ind1, ind2) in itertools.combinations(index_list, 2)\
								if (Z.index(ind, ind).level()[0]) > 0.5]

		graph 							= nx.Graph()

		graph.add_nodes_from(self.vertices)
		graph.add_edges_from(z_star_edges)

		components 						= [list(component) for component in nx.connected_components(graph)]

		unassigned_partitions 			= deepcopy(self.partitions)
		unassigned_partitions.remove(self.graph.nodes[self.fixed_vertex]["partition"])

		for (partition, component) in enumerate(components):

			if self.fixed_vertex in component:
				for vertex in component:
					self.graph.nodes[vertex]["partition"] 	= self.graph.nodes[self.fixed_vertex]["partition"]
			else:
				for vertex in component:
					self.graph.nodes[vertex]["partition"] 	= unassigned_partitions[0]
				unassigned_partitions.remove(unassigned_partitions[0])


	else: 
		sys.exit("Mosek solver cannot find the coefficients in the specified time-limit.")



	self.gurobi_end_time	 		= time.time()
	self.print_gurobi_results_summary()
					
	
	
	






