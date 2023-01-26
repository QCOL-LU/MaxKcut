import networkx as nx 
from gurobipy import *
import time
import itertools
from copy import deepcopy




#================================================================================================
# The Partition-based MILO formulation
#================================================================================================
def solve_max_k_cut_pmilo(self):

	self.gurobi_start_time 			= time.time()

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

			z 			= {(vertex1, vertex2): model.addVar(vtype=GRB.BINARY, name="z(%i,%i)" %(vertex1, vertex2)) \
							for (vertex1, vertex2) in itertools.combinations(self.vertices, 2)}

			model._z 	= z

			#-----------------------------------------------------------------------------------
			# Objective function
			#-----------------------------------------------------------------------------------
			model.setObjective(quicksum(self.graph.edges[edge]["weight"] * (1 - z[edge]) for edge in self.edges), GRB.MAXIMIZE)
			

			#-----------------------------------------------------------------------------------
			# Constraints
			#-----------------------------------------------------------------------------------
			for vertex_set in itertools.combinations(self.vertices, 3):
				model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[0], vertex_set[2])] <= 1 + z[(vertex_set[1], vertex_set[2])])
				model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[2])])
				model.addConstr(z[(vertex_set[0], vertex_set[2])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[1])])


			if self.Params.Relaxed == True:
				for clique in itertools.combinations(self.vertices, self.num_partitions + 1):
					clique_list 			= sorted(list(clique))

					model.addConstr(quicksum(z[edge] for edge in itertools.combinations(clique_list, 2) ) >=  1 )
					# clique_constraint.Lazy 	= -1

			#-----------------------------------------------------------------------------------
			# Callback function to obtain root node relaxation
			#-----------------------------------------------------------------------------------
			def add_lazy_constraint(model, where):
				if where == GRB.Callback.MIPSOL:

					z_star_edges 		= [(vertex1, vertex2) for (vertex1, vertex2) in itertools.combinations(self.vertices, 2)\
											 if (model.cbGetSolution(model._z[(vertex1, vertex2)]) ) > .5]

					graph 				= nx.Graph()

					graph.add_nodes_from(self.vertices)
					graph.add_edges_from(z_star_edges)

					components 			= [list(component)[0] for component in nx.connected_components(graph)]

					size 				= len(components)

					if size <= self.num_partitions:
						return


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
				model.optimize()
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
			z_star_edges 					= [(vertex1, vertex2) for (vertex1, vertex2) in itertools.combinations(self.vertices, 2)\
										 		if (z[(vertex1, vertex2)].x) > 0.5]

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


	self.gurobi_end_time	 		= time.time()
	self.print_gurobi_results_summary()
					
	

