import networkx as nx 
from gurobipy import *
import time



#================================================================================================
# Solve the BQO formulation
#================================================================================================
def solve_max_k_cut_qubo(self):

	self.gurobi_start_time 			= time.time()

	penalty_coef	= {vertex: max(self.graph.nodes[vertex]["pos-weight"]/self.num_partitions, -self.graph.nodes[vertex]["neg-weight"]/2) for vertex in self.vertices}
	
	# penalty_coef[1] = 1
	penalty_coef[1] -= 1e-2
	print(penalty_coef)

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

			self.fixed_vertex = 3
			self.graph.nodes[self.fixed_vertex]["partition"] = 0


			#-----------------------------------------------------------------------------------
			# Objective function
			#-----------------------------------------------------------------------------------

			model.setObjective(	self.total_weights - quicksum(self.graph.edges[edge]["weight"] * x[edge[0],partition]*x[edge[1],partition] for edge in self.edges for partition in self.partitions) \
								- quicksum((penalty_coef[vertex]) * (quicksum(x[vertex, partition] for partition in self.partitions) -  1) * (quicksum(x[vertex, partition] for partition in self.partitions) -  1) for vertex in self.vertices), GRB.MAXIMIZE)		
			
			

			model.update()
			

			for partition in self.partitions:
				if partition == self.graph.nodes[self.fixed_vertex]["partition"]:
					x[self.fixed_vertex, partition].lb 	= 1
				else:
					x[self.fixed_vertex, partition].ub 	= 0


			model.update()
			

			#-----------------------------------------------------------------------------------
			# Set the solver parameters
			#-----------------------------------------------------------------------------------
			model.update()
			model.Params.Heuristics 		= self.Params.Gurobi_Heuristics
			model.Params.Presolve 			= self.Params.Gurobi_Presolve
			model.Params.Symmetry 			= self.Params.Gurobi_Symmetry
			model.Params.Cuts 				= self.Params.Gurobi_Cuts

			# model.Params.NonConvex 			= -1 if self.Params.Curvature_Type == "convex" else 2
			model.Params.NonConvex 			= 2
			
			model.Params.Threads			= self.Params.Gurobi_Threads
			model.Params.timeLimit 			= self.Params.Gurobi_TimeLimit

			# model.Params.LogToConsole 		= 0
			model.Params.LogFile 			= self.filename
			# model.Params.OutputFlag 		= 0
			
			model 							= model if self.Params.Relaxed == False else model.relax()
			model.update()

			
			#-----------------------------------------------------------------------------------
			# Solve the model and extract the obtained solution
			#-----------------------------------------------------------------------------------
			model.optimize()

			self.gurobi_obj_value			= model.objVal
			self.gurobi_BB_nodes 			= model.getAttr(GRB.Attr.NodeCount) 		# Spatial BB nodes
			self.gurobi_ObjBound			= model.ObjBound
			self.gurobi_MIPGap				= model.MIPGap

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
			print(solution)
			print(hkj)
			
			vertices_zero_assignment 					= [vertex for vertex in self.vertices if sum(solution[vertex]) == 0]
			vertices_more_than_one_assignment			= [vertex for vertex in self.vertices if sum(solution[vertex]) > 1]
			
			for vertex in vertices_zero_assignment:
				num_assigned_partition 					= sum(solution[vertex])
				weight_partition_based_neighbors 		= {partition: sum([self.graph.edges[vertex, neighbor]["weight"] for neighbor in self.graph.neighbors(vertex) if solution[neighbor][partition] == 1]) for partition in self.partitions}
				selected_partition 						= min(weight_partition_based_neighbors, key=weight_partition_based_neighbors.get)
				solution[vertex] 						= [1 if partition == selected_partition else 0 for partition in self.partitions]

				print(vertex, "z")
			while vertices_more_than_one_assignment:
				vertex 	 								= vertices_more_than_one_assignment[0]
				similar_inf_vertices					= [neighbor for neighbor in vertices_more_than_one_assignment[1:] if str(solution[vertex]) == str(solution[neighbor]) ]
				inf_edges 								= [edge for edge in self.edges if self.graph.edges[edge]["weight"] < 0 and len(set(edge).intersection(set(similar_inf_vertices))) > 1 ]
				inf_partitions 							= [partition for partition in self.partitions if solution[vertex][partition] == 1]
				weight_partition 						= {partition: sum(self.graph.edges[edge[0], edge[1]]["weight"] * solution[edge[0]][partition] * solution[edge[1]][partition] for edge in inf_edges) for partition in self.partitions}

				selected_partition 						= min(weight_partition, key=weight_partition.get)
				solution[vertex] 						= [1 if partition == selected_partition else 0 for partition in self.partitions]
				vertices_more_than_one_assignment 		= list(set(vertices_more_than_one_assignment) - set(similar_inf_vertices))
				print(vertex, "o")

	for vertex in self.vertices:
		for partition in self.partitions:

			if solution[vertex][partition] > .5:
				self.graph.nodes[vertex]["partition"] 	= partition 

	self.gurobi_end_time	 		= time.time()

	self.print_gurobi_results_summary()
					
	

	