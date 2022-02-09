import networkx as nx 
import mosek.fusion as mosek
import numpy as np
import scipy as sp
import time
import threading
import itertools
import sys



#================================================================================================
# Calculate coefficients to convexify/concavify the BQO model
#================================================================================================
def calculate_curvature_coefs(self):

	self.start_time_curv				= time.time()
	adjacency_matrix 					= nx.adjacency_matrix(self.graph, weight='weight').todense()
	
	if self.Params.Method == "BQO":

		self.variable_coefs 			= {vertex:0 for vertex in self.vertices}

		hessian_matrix_dim 				= self.num_vertices
		hessian_matrix 					= adjacency_matrix
		

	elif self.Params.Method == "R-BQO":

		self.variable_coefs  			= {vertex:{partition:0  for partition in self.partitions[:-1]} for vertex in self.vertices}
		index_convertor 				= {partition * self.num_vertices + ind: (vertex, partition) for partition in self.partitions[:-1] for (ind, vertex) in enumerate(self.vertices)}

		hessian_matrix_dim 				= self.num_vertices * (self.num_partitions - 1)
		hessian_matrix 					= np.zeros((hessian_matrix_dim, hessian_matrix_dim))
		
		for (ind1, ind2) in itertools.combinations(range(hessian_matrix_dim), 2):

			(vertex1, partition1)		= index_convertor[ind1]
			(vertex2, partition2)		= index_convertor[ind2]

			adj_ind1					= int(ind1/self.num_partitions)
			adj_ind2					= int(ind2/self.num_partitions)
			
			hessian_matrix[ind1, ind2] 	= adjacency_matrix[adj_ind1, adj_ind2] * (2 if partition1 == partition2 else 1)
			hessian_matrix[ind2, ind1] 	= adjacency_matrix[adj_ind1, adj_ind2] * (2 if partition1 == partition2 else 1)



		
	#--------------------------------------------------------------------------------------------
	# Using MOSEK solver to determine the exact eigenvalues of the adjacency matrix
	#--------------------------------------------------------------------------------------------
	if self.Params.Curvature_Method == "mosek":
		#----------------------------------------------------------------------------------------
		# Model initialization
		#----------------------------------------------------------------------------------------
		model 			= mosek.Model("curvature coefficients")
		D_plus_A 		= model.variable(mosek.Domain.isTrilPSD(hessian_matrix_dim))

		
		#----------------------------------------------------------------------------------------
		# Objective function
		#----------------------------------------------------------------------------------------
		objective 		= D_plus_A.index(0, 0)
		for ind in range(1, hessian_matrix_dim):
			objective 	= mosek.Expr.add(objective, D_plus_A.index(ind, ind))

		model.objective(mosek.ObjectiveSense.Minimize, objective)

	
		#----------------------------------------------------------------------------------------
		# Constraints
		#----------------------------------------------------------------------------------------
		for (col, row) in itertools.combinations(range(hessian_matrix_dim), 2):
			model.constraint(D_plus_A.index(col, row), mosek.Domain.equalsTo(adjacency_matrix[col, row]))
			model.constraint(D_plus_A.index(row, col), mosek.Domain.equalsTo(adjacency_matrix[row, col]))


		#----------------------------------------------------------------------------------------
		# Solve the SDO model
		#----------------------------------------------------------------------------------------
		time_limit 				= 60
		
		T 		= threading.Thread(target=model.solve)

		try:
			T.start() 			# optimization now running in background

			# Loop until we get a solution or you run out of patience and press Ctrl-C
			while True:
				if not T.is_alive():
					print("Solver terminated before anything happened!")
					break
				elif time.time() - self.start_time_curv > time_limit:
					print("Solver terminated due to time_limit!")
					model.breakSolver()
					break
		except KeyboardInterrupt:
			print("Signalling the solver that it can give up now!")
			model.breakSolver()
		finally:
			try: T.join() 		# wait for the solver to return
			except: pass


		#----------------------------------------------------------------------------------------
		# Extract the obtained solution
		#----------------------------------------------------------------------------------------
		if model.getPrimalSolutionStatus() == mosek.SolutionStatus.Feasible or model.getPrimalSolutionStatus() == mosek.SolutionStatus.Optimal:
			
			if self.Params.Method == "BQO":

				for (ind, variable) in enumerate(self.variable_coefs.keys()):
					self.variable_coefs[variable] 			= D_plus_A.index(ind, ind).level()[0] * (-1 if self.Params.Curvature_Type == "concave" else 1)
			

			elif self.Params.Method == "R-BQO":

				for ind in enumerate(range(hessian_matrix_dim)):
					(vertex, partition)						= index_convertor(ind)
					self.variable_coefs[vertex][partition]	= D_plus_A.index(ind, ind).level()[0] * (-1 if self.Params.Curvature_Type == "concave" else 1)

		else: 
			sys.exit("Mosek solver cannot find the coefficients in the specified time-limit.")



	#--------------------------------------------------------------------------------------------
	# Using Power Iteration Method to determine the approximate eigenvalues of the adjacency matrix
	#--------------------------------------------------------------------------------------------
	elif self.Params.Curvature_Method == "power_iteration": 
		
		max_eigenvalue_vector 				= np.matrix([np.ones(hessian_matrix_dim)]).T
		max_eigenvalue_vector_old			= 2 * np.matrix([np.ones(hessian_matrix_dim)]).T

		if self.Params.Method == "BQO":
			eigenvalue_lower_bound 			= - adjacency_matrix.max() * np.sqrt(hessian_matrix_dim/2 * np.floor((hessian_matrix_dim + 1)/2) )
			eigenvalue_upper_bound  		= - eigenvalue_lower_bound * (hessian_matrix_dim - 1)
		
		elif self.Params.Method == "R-BQO": 
			eigenvalue_upper_bound 			= max(hessian_matrix.sum(axis=0))		# Based on Perron-Frobenius Theorem
			eigenvalue_lower_bound 			= - eigenvalue_upper_bound


		if self.Params.Curvature_Type == "convex":
			modified_matrix					= eigenvalue_upper_bound * np.eye(hessian_matrix_dim) - hessian_matrix

		elif self.Params.Curvature_Type == "concave":
			modified_matrix					= hessian_matrix - eigenvalue_lower_bound * np.eye(hessian_matrix_dim) 
		

		while np.linalg.norm(max_eigenvalue_vector - max_eigenvalue_vector_old) >  1e-3:
			max_eigenvalue_vector_old 		= max_eigenvalue_vector
			max_eigenvalue_vector 			= max_eigenvalue_vector / np.linalg.norm(max_eigenvalue_vector)
			max_eigenvalue_vector 			= np.matmul(modified_matrix, max_eigenvalue_vector)


		if self.Params.Curvature_Type == "convex":
			min_eigenvalue 					= eigenvalue_upper_bound - max(max_eigenvalue_vector)	

		elif self.Params.Curvature_Type == "concave":
			max_eigenvalue 					= eigenvalue_lower_bound + max(max_eigenvalue_vector)	


		if self.Params.Method == "BQO":
			for variable in self.variable_coefs.keys():
				self.variable_coefs[variable] 				= - (max_eigenvalue if self.Params.Curvature_Type == "concave" else min_eigenvalue)

		elif self.Params.Method == "R-BQO":
			for ind in range(hessian_matrix_dim):
				(vertex, partition)							= index_convertor[ind]
				self.variable_coefs[vertex][partition]		= - (max_eigenvalue if self.Params.Curvature_Type == "concave" else min_eigenvalue)



	#--------------------------------------------------------------------------------------------
	# Using NUMPY Python package to determine the exact eigenvalues of the adjacency matrix
	#--------------------------------------------------------------------------------------------
	elif self.Params.Curvature_Method == "numpy":
		eigenvalues 						= np.linalg.eigvals(hessian_matrix)

		if self.Params.Method == "BQO":
			for (ind, variable) in enumerate(self.variable_coefs.keys()):
				self.variable_coefs[variable] 				= eigenvalues[ind].real * (1 if self.Params.Curvature_Type == "concave" else -1)

		elif self.Params.Method == "R-BQO":
			for ind in range(hessian_matrix_dim):
				(vertex, partition)							= index_convertor[ind]
				self.variable_coefs[vertex][partition]		= eigenvalues[ind].real * (1 if self.Params.Curvature_Type == "concave" else -1)

	elif self.Params.Curvature_Method == "bound":
		if self.Params.Method == "BQO":
			for (ind, variable) in enumerate(self.variable_coefs.keys()):
				self.variable_coefs[variable] 				= np.sqrt(self.num_vertices/2 * np.floor( (self.num_vertices + 1)/2 ))

	self.end_time_curv 						= time.time()
