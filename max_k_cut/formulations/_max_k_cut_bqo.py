import networkx as nx 
from gurobipy import *
import time
import itertools
from copy import deepcopy



#================================================================================================
# Solve the BQO formulation
#================================================================================================
def solve_max_k_cut_bqo(self):

	self.gurobi_start_time 			= time.time()


	if self.Params.Clique_Constraints == True:
		#--------------------------------------------------------------------------------------------
		# Construct the extended chordal graph
		#--------------------------------------------------------------------------------------------
		chordal_graph, 	H					= nx.complete_to_chordal_graph(self.graph)
		maximal_cliques 					= list(nx.find_cliques(chordal_graph))
		chordal_graph_edges 				= [ (min(edge), max(edge)) for edge in chordal_graph.edges()]
		
		self.density_chordal_graph 			= 200 * chordal_graph.number_of_edges() / ( (self.num_vertices - 1) * self.num_vertices )


	#--------------------------------------------------------------------------------------------
	# Model initialization
	#--------------------------------------------------------------------------------------------
	with Env(empty=True) as env:
		if self.Params.Gurobi_LogToConsole == 0:
			env.setParam('LogToConsole', 0)
			sys.stdout.write(2*"\033[F\033[K")
		env.start()
		with Model(env=env) as model:			

			#-----------------------------------------------------------------------------------
			# Variable definition
			#-----------------------------------------------------------------------------------
			x 		= model.addVars(self.vertices, self.partitions, vtype=GRB.BINARY, name="x")
			model._x 	= x

			# model.update()

			if self.Params.Symmetry_Breaking == True:
				s 	= model.addVars(self.vertices, self.partitions, vtype=GRB.CONTINUOUS) 
				u 	= model.addVars(self.vertices, self.partitions, vtype=GRB.CONTINUOUS) 
				w 	= model.addVars(self.vertices, self.partitions, vtype=GRB.CONTINUOUS)
				r 	= model.addVars(self.vertices, self.partitions, vtype=GRB.BINARY)

			
			#-----------------------------------------------------------------------------------
			# Objective function
			#-----------------------------------------------------------------------------------	
			if self.Params.Curvature_Type == "indefinite":

				model.setObjective(	self.total_weights - quicksum(self.graph.edges[edge]["weight"] * x[edge[0],partition]*x[edge[1],partition] for edge in self.edges for partition in self.partitions) \
									, GRB.MAXIMIZE)		
			
			elif self.Params.Curvature_Type in ["convex", "concave"]:

				self.calculate_curvature_coefs()

				model.setObjective(	self.total_weights - quicksum(self.graph.edges[edge]["weight"] * x[edge[0],partition]*x[edge[1],partition] for edge in self.edges for partition in self.partitions) \
									+ quicksum(self.variable_coefs[vertex] * (x[vertex, partition] *  x[vertex, partition] - x[vertex, partition]) \
									 for vertex in self.vertices for partition in self.partitions), GRB.MAXIMIZE)
			else:
				sys.exit("Please specify the correct objective type: (indefinite, convex, concave)")


			model.update()
			
			#-----------------------------------------------------------------------------------
			# Constraints
			#-----------------------------------------------------------------------------------
			for vertex in self.vertices:
				model.addConstr(quicksum(x[vertex, partition] for partition in self.partitions) == 1)

				if self.Params.Relaxed == True and self.Params.Relaxed_NonConvex == True:
					for partition in self.partitions:
						model.addConstr(x[vertex, partition] <= x[vertex, partition]*x[vertex, partition])



			for partition in self.partitions:
				if partition == self.graph.nodes[self.fixed_vertex]["partition"]:
					x[self.fixed_vertex, partition].lb 	= 1
				else:
					x[self.fixed_vertex, partition].ub 	= 0


			if self.Params.Clique_Constraints == True:
				for clique in maximal_cliques:
					clique_size 	= len(clique)
					(quotient, reminder)	= divmod(clique_size, self.num_partitions)
					# if clique_size < self.num_partitions + 1: continue #
					model.addConstr(quicksum( (x[vertex1, partition] - x[vertex2, partition]) * (x[vertex1, partition] - x[vertex2, partition])
												for (vertex1, vertex2) in itertools.combinations(clique, 2)
												for partition in self.partitions) 
									<= ((clique_size*(clique_size - quotient)) - reminder * (quotient + 1) )  )


			model.update()
			#-----------------------------------------------------------------------------------
			# Symmetry Breaking Constraints
			#-----------------------------------------------------------------------------------
			if self.Params.Symmetry_Breaking == True:
				# line 47
				r[self.vertices[0], 0].lb 	= 1

				for (ind, vertex) in enumerate(self.vertices):

					for partition in self.partitions:
						
						if partition == self.num_partitions - 1:
							# line 35
							model.addConstr(x[vertex, partition] == s[vertex, partition])

							if ind == len(self.vertices) - 1:
								# line 48
								model.addConstr(u[vertex, partition] + r[vertex, partition] == 1 )

							else:
								# line 44
								model.addConstr(u[vertex, partition] + r[vertex, partition] == u[self.vertices[ind + 1], partition] )
								
						else:
							# line 34
							model.addConstr(x[vertex, partition] == s[vertex, partition] - s[vertex, partition + 1])

							if ind == len(self.vertices) - 1:
								# line 45
								model.addConstr(u[vertex, partition] + r[vertex, partition] == 0 )
							else:
								# line 43
								model.addConstr(u[vertex, partition] + r[vertex, partition] == u[self.vertices[ind + 1], partition] + r[self.vertices[ind + 1], partition + 1] )			


						# lines 37 38
						if ind == 0:
							model.addConstr(r[vertex, partition] == w[vertex, partition])
						else:
							model.addConstr(r[vertex, partition] == w[vertex, partition] - w[self.vertices[ind - 1], partition])


						# lines 40 41
						model.addConstr(r[vertex, partition] <= x[vertex, partition])
						model.addConstr(s[vertex, partition] <= w[vertex, partition])


			#-----------------------------------------------------------------------------------
			# Callback function to obtain root node relaxation
			#-----------------------------------------------------------------------------------
			def add_lazy_constraint(model, where):
				if where == GRB.Callback.MIPSOL:

					for clique in maximal_cliques:
						if len(clique) <= self.num_partitions:
							continue
							
						clique.sort()

	
						z_star_edges 		= [(vertex1, vertex2) for (vertex1, vertex2) in itertools.combinations(clique, 2) \
												 if sum(model.cbGetSolution(model._x[vertex1, partition]) 
												 		* model.cbGetSolution(model._x[vertex2, partition]) 
												 		for partition in self.partitions) > .5]

						graph 				= nx.Graph()

						graph.add_nodes_from(clique)
						graph.add_edges_from(z_star_edges)

						components 			= [list(component)[0] for component in nx.connected_components(graph)]

						size 				= len(components)

						if size <= self.num_partitions:
							continue


						for (ind, _) in enumerate(components):
							first 			= components[ind: min(size, ind + self.num_partitions + 1)]
							size_second 	= self.num_partitions + 1 - len(first)
							second 			= components[:size_second]
							set_Q 			= first + second

							model.cbLazy(quicksum((model._x[key1, partition] - model._x[key2, partition])* (model._x[key1, partition] - model._x[key2, partition])
												for (ind1, key1) in enumerate(set_Q[:-1]) 
												for key2 in set_Q[ind1 + 1:] 
												for partition in self.partitions) 
										<= (self.num_partitions + 2) * (self.num_partitions + 1) - 2)
							model.update()

			#-----------------------------------------------------------------------------------
			# Set the solver parameters
			#-----------------------------------------------------------------------------------
			# model.update()
			model.Params.Heuristics 		= self.Params.Gurobi_Heuristics
			model.Params.Presolve 			= self.Params.Gurobi_Presolve
			model.Params.Symmetry 			= self.Params.Gurobi_Symmetry
			model.Params.Cuts 				= self.Params.Gurobi_Cuts
			# model.Params.RLTCuts 			= 2

			model.Params.NonConvex 			= -1 if self.Params.Curvature_Type == "concave" else 2
			
			model.Params.Threads			= self.Params.Gurobi_Threads
			model.Params.timeLimit 			= self.Params.Gurobi_TimeLimit - ( time.time() - self.gurobi_start_time)
			model.Params.MIPGap 			= self.Params.Gurobi_MIPGap
			# model.Params.LazyConstraints 	= 1 if self.Params.Clique_Constraints else 0


			# model.Params.LogToConsole 		= 0
			model.Params.LogFile 			= self.filename[:-4] + "_log.txt"
			# model.Params.OutputFlag 		= 0
			
			if self.Params.Relaxed == True:
				model 						= model.relax()
			# model.update()

			
			#-----------------------------------------------------------------------------------
			# Solve the model and extract the obtained solution
			#-----------------------------------------------------------------------------------
			# if self.Params.Clique_Constraints == True:
			# 	model.optimize(add_lazy_constraint)
			# else:
			model.optimize()

			self.gurobi_obj_value			= model.objVal
			self.gurobi_BB_nodes 			= model.getAttr(GRB.Attr.NodeCount) 		# Spatial BB nodes
			if self.Params.Relaxed == False or self.Params.Curvature_Type != "concave":
				self.gurobi_ObjBound			= model.ObjBound
				self.gurobi_MIPGap				= model.MIPGap
			else:
				self.gurobi_ObjBound			= model.objVal
				self.gurobi_MIPGap				= 0.0

			

			if model.status == GRB.Status.OPTIMAL:
				self.gurobi_model_status 	= "optimal"
			
			elif model.status == GRB.Status.INFEASIBLE:
				self.gurobi_model_status 	= "infeasible"

			elif model.status == GRB.Status.TIME_LIMIT:
				self.gurobi_model_status 	= "time limit"

			elif model.status == GRB.Status.INTERRUPTED:
				self.gurobi_model_status 	= "interrupted"

			else:
				self.gurobi_model_status 	= model.status 
		

			#-----------------------------------------------------------------------------------
			# Obtain a binary solution
			#-----------------------------------------------------------------------------------
			solution						= {vertex: [model.getVarByName("x["+str(vertex)+","+str(partition)+"]").x for partition in self.partitions] for vertex in self.vertices}

	if self.Params.Rounding_Heuristic == True and self.Params.Relaxed == True: 
		for vertex in self.vertices:
			b_vector 				= [sum([self.graph.edges[(vertex, neighbor)]["weight"] * solution[neighbor][partition] for neighbor in self.graph.neighbors(vertex)]) for partition in self.partitions]
			selected_partition 		= b_vector.index(min(b_vector))
			solution[vertex] 		= [1 if partition == selected_partition else 0 for partition in self.partitions]
			

	if self.Params.Rounding_Heuristic == True or self.Params.Relaxed == False:
		
		for vertex in self.vertices:
			for partition in self.partitions:

				if solution[vertex][partition] > .5:
					self.graph.nodes[vertex]["partition"] 	= partition 

	self.gurobi_end_time	 		= time.time()

	self.print_gurobi_results_summary()
					
	

	