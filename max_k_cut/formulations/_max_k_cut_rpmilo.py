import networkx as nx 
from gurobipy import *
import time
import itertools
from copy import deepcopy
from random import shuffle
import os, psutil


from networkx.algorithms.chordal import complete_to_chordal_graph


#================================================================================================
# The Reduced Partition-based MILO formulation proposed by Wang and Hijazi
#================================================================================================
def solve_max_k_cut_rpmilo(self):

	self.gurobi_start_time 			= time.time()

	#--------------------------------------------------------------------------------------------
	# Construct the extended chordal graph
	#--------------------------------------------------------------------------------------------
	chordal_graph, 	H					= complete_to_chordal_graph(self.graph)
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
			process = psutil.Process(os.getpid())
	
			#-----------------------------------------------------------------------------------
			# Variable definition
			#-----------------------------------------------------------------------------------
			# z 			= {(vertex1, vertex2): model.addVar(vtype=GRB.BINARY, name="z(%i,%i)" %(vertex1, vertex2)) \
			# 				for (vertex1, vertex2) in chordal_graph_edges}

			z 			= model.addVars(chordal_graph_edges, vtype=GRB.BINARY, name="z")

			model._z 	= z

			#-----------------------------------------------------------------------------------
			# Objective function
			#-----------------------------------------------------------------------------------
			model.setObjective(quicksum(self.graph.edges[edge]["weight"] * (1 - z[edge]) for edge in self.edges), GRB.MAXIMIZE)
			

			#-----------------------------------------------------------------------------------
			# Constraints
			#-----------------------------------------------------------------------------------
			for clique in maximal_cliques:
				len_clique			= len(clique)
				if len_clique >= 3:
					for vertex_set in itertools.combinations(clique, 3):
						vertex_set 			= sorted(vertex_set)
						model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[0], vertex_set[2])] <= 1 + z[(vertex_set[1], vertex_set[2])])
						model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[2])])
						model.addConstr(z[(vertex_set[0], vertex_set[2])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[1])])

				# diff			= len_clique - (self.num_partitions + 1)
				# if diff >= 0:
				# 	number 		= (len_clique - 1) if diff >= 1 else (self.num_partitions + 1)
					
				# 	for vertex_set in itertools.combinations(clique, number):
				# 		vertex_set 			= sorted(vertex_set)
				# 		model.addConstr(quicksum(z[edge] for edge in itertools.combinations(vertex_set, 2) ) >=  1 )

				# if self.Params.Relaxed == True:
				# 	for subclique in itertools.combinations(clique, self.num_partitions + 1):
				# 		clique_list 			= sorted(list(subclique))

				# 		clique_constraint 		= model.addConstr(quicksum(z[edge] for edge in itertools.combinations(clique_list, 2) ) >=  1 )
				# 		clique_constraint.Lazy 	= 3

				# if len(clique) >= self.num_partitions + 1:
				# 	for vertex_set in itertools.combinations(clique, self.num_partitions + 1):
		
				# 		clique_constraint		= model.addConstr(quicksum(z[(min(vertex1, vertex2), max(vertex1, vertex2) )] for (vertex1, vertex2) in itertools.combinations(vertex_set, 2) ) >= 1 )
				# 		clique_constraint.Lazy 	= 3

			# if self.Params.Clique_Constraints == True:
			# 	for clique in nx.find_cliques(self.graph):
			# 		clique_list 			= sorted(list(clique))
			# 		clique_size 			= len(clique_list)
			# 		(quotient, reminder)	= divmod(clique_size, self.num_partitions)

			# 		clique_constraint 		= model.addConstr(quicksum(z[edge] for edge in itertools.combinations(clique_list, 2) ) >=  ((quotient*(quotient - 1)*(self.num_partitions - reminder) + quotient*(quotient + 1)*reminder) )/2 )
			# 		clique_constraint.Lazy 	= -1

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
												 if (model.cbGetSolution(model._z[(vertex1, vertex2)]) ) > .5]

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

							model.cbLazy(quicksum(model._z[(min(key1, key2), max(key1, key2) )] \
									for (ind1, key1) in enumerate(set_Q[:-1]) for key2 in set_Q[ind1 + 1:] ) >= 1)
							model.update()


				

			#-----------------------------------------------------------------------------------
			# Set the solver parameters
			#-----------------------------------------------------------------------------------
			model.update()
			model.Params.Heuristics 		= self.Params.Gurobi_Heuristics
			model.Params.Presolve 			= self.Params.Gurobi_Presolve
			model.Params.Symmetry 			= self.Params.Gurobi_Symmetry
			model.Params.Cuts 				= self.Params.Gurobi_Cuts
			
			model.Params.Threads			= self.Params.Gurobi_Threads
			model.Params.timeLimit 			= self.Params.Gurobi_TimeLimit
			model.Params.MIPGap 			= self.Params.Gurobi_MIPGap 

			model.Params.LogFile 			= self.filename[:-4] + "_log.txt"
			model.Params.LazyConstraints 	= 1

			model 							= model.relax() if self.Params.Relaxed == True else model
			model.update()

			#-----------------------------------------------------------------------------------
			# Solve the model and extract the obtained solution
			#-----------------------------------------------------------------------------------
		
			if self.Params.Relaxed == True:
				model.Params.LPWarmStart 			= 1

				model.optimize()
				used_memory_gib 		= process.memory_info().rss / 1024 ** 3
				print("Used memory: {:.02f}GB".format(used_memory_gib) )
				
				available_memory 		= psutil.virtual_memory()[1] >> 30
				limit 					= 100000 // len(maximal_cliques)
				if available_memory > 1: 
					done 		= False
					threshold 	= 0.1
					max_rhs 	= (1 - 1e-6)

					def reached_time_limit():
						return time.time() - self.gurobi_start_time > self.Params.Gurobi_TimeLimit 


					while (not done) and (available_memory > 5) and (not reached_time_limit()):
						# print("threshold:", threshold)

						done 		= True
						num_constrs		= 1
					
						for edge in chordal_graph_edges:
							chordal_graph.edges[edge]["weight"] = model.getVarByName(z[edge].VarName).x

						for clique in maximal_cliques:

							# len_clique				= len(clique)
							
							# if len_clique >= 3:
							# 	num_constrs				= 1
								
							# 	for vertex_set in itertools.combinations(clique, 3):
							# 		vertex_set 			= sorted(vertex_set)
							# 		edge0, edge1, edge2 = list(itertools.combinations(vertex_set, 2))

							# 		not_satsified0 		= chordal_graph.edges[edge2]["weight"] + chordal_graph.edges[edge1]["weight"] - (chordal_graph.edges[edge0]["weight"] + 1) < threshold 
							# 		not_satsified1 		= chordal_graph.edges[edge0]["weight"] + chordal_graph.edges[edge2]["weight"] - (chordal_graph.edges[edge1]["weight"] + 1) < threshold 
							# 		not_satsified2 		= chordal_graph.edges[edge0]["weight"] + chordal_graph.edges[edge1]["weight"] - (chordal_graph.edges[edge2]["weight"] + 1) < threshold 
								
							# 		if reached_time_limit(): break
							# 		if not_satsified0:
							# 			model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[0], vertex_set[2])] <= 1 + z[(vertex_set[1], vertex_set[2])])
							# 			num_constrs 		+= 1

							# 		if not_satsified1:	
							# 			model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[2])])
							# 			num_constrs 		+= 1

							# 		if not_satsified2:
							# 			model.addConstr(z[(vertex_set[0], vertex_set[2])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[1])])
							# 			num_constrs 		+= 1

							# 		if num_constrs > limit: break

							if len(clique) >= self.num_partitions + 1:

								subgraph 		= chordal_graph.subgraph(clique)
								clique_weight 	= {vertex: subgraph.degree(vertex, weight="weight") for vertex in clique}

								sorted_clique 	= sorted(clique_weight, key=lambda k: clique_weight[k])

								num_constrs		= 1

								for subclique in itertools.combinations(sorted_clique, self.num_partitions + 1):

									sorted_subclique 		= sorted(list(subclique))
									not_satsified 			= sum(chordal_graph.edges[edge]["weight"] for edge in itertools.combinations(sorted_subclique, 2) ) < threshold 
									
									if reached_time_limit(): break
									if not_satsified:

										model.addConstr(quicksum(model.getVarByName(z[edge].VarName) for edge in itertools.combinations(sorted_subclique, 2) ) >=  1 )

										done 				= False

										if num_constrs > limit: break
										
										num_constrs 		+= 1

							if num_constrs > limit: break
						
						available_memory 		= psutil.virtual_memory()[1] >> 30
						if reached_time_limit(): break	
						
						if not done: 
							model.Params.timeLimit = self.Params.Gurobi_TimeLimit - (time.time() - self.gurobi_start_time)
							model.optimize()
						
						elif threshold < max_rhs:
							threshold 				= min(max_rhs, threshold + 0.1)
							done 					= False


							
						

			else:
				model.optimize(add_lazy_constraint)

			self.gurobi_obj_value			= model.objVal
			self.gurobi_BB_nodes 			= model.getAttr(GRB.Attr.NodeCount) 		# Spatial BB nodes
			self.gurobi_ObjBound			= model.ObjBound
			self.gurobi_MIPGap				= model.MIPGap if self.Params.Relaxed == False else 0.0
			
			all_variables 					= model.getVars()

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
			if self.Params.Relaxed == False:
				z_star_edges 				= [(vertex1, vertex2) for (vertex1, vertex2) in chordal_graph_edges if (z[(vertex1, vertex2)].x) > 0.5]

				conflict_vertex_pairs 		= [(vertex1, vertex2) for (vertex1, vertex2) in chordal_graph_edges if (z[(vertex1, vertex2)].x) <= 0.5]

				graph 							= nx.Graph()

				graph.add_nodes_from(self.vertices)
				graph.add_edges_from(z_star_edges)

				components 						= [list(component) for component in nx.connected_components(graph)]
				num_components 					= len(components)

				vertex_components 				= {vertex: -1 for vertex in self.vertices}
				
				for vertex in self.vertices: 
					for (ind, component) in enumerate(components):
						if vertex in component:
							vertex_components[vertex] 	=  ind
							break

				with Model(env=env) as model:
					x 		= model.addVars(range(num_components), self.partitions, vtype=GRB.BINARY, name="x")

					model.addConstrs(quicksum(x[ind, partition] for partition in self.partitions) == 1 for ind in range(num_components)) 

					x[vertex_components[self.fixed_vertex], self.graph.nodes[self.fixed_vertex]["partition"]].lb 		= 1

					for (vertex1, vertex2) in conflict_vertex_pairs:
						ind1  	= vertex_components[vertex1]
						ind2  	= vertex_components[vertex2]

						model.addConstrs(x[ind1, partition] + x[ind2, partition] <= 1 for partition in self.partitions) 


					model.optimize()

					for (ind, component) in enumerate(components):
						partition  		= [partition for partition in self.partitions if x[ind, partition].x > .5][0]
						for vertex in component:
							self.graph.nodes[vertex]["partition"] 	= partition


	self.gurobi_end_time	 		= time.time()
	self.print_gurobi_results_summary()
					
	

