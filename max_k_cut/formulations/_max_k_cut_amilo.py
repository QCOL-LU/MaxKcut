import networkx as nx 
from gurobipy import *
import time
import itertools
import numpy as np


#================================================================================================
# The Assignment-based MILO formulation
#================================================================================================
def solve_max_k_cut_amilo(self):

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
			x 			= model.addVars(self.vertices, self.partitions, vtype=GRB.BINARY, name="x")
			y 			= model.addVars(self.edges, vtype=GRB.BINARY)

			#-----------------------------------------------------------------------------------
			# Objective function
			#-----------------------------------------------------------------------------------
			model.setObjective( quicksum(self.graph.edges[edge]["weight"] * y[edge] for edge in self.edges), GRB.MAXIMIZE)
			
			#-----------------------------------------------------------------------------------
			# Constraints
			#-----------------------------------------------------------------------------------
			for vertex in self.vertices:
				model.addConstr(quicksum(x[vertex, partition] for partition in self.partitions) == 1)
				
			for partition in self.partitions:
				for edge in self.edges:

					model.addConstr( x[edge[0], partition] - x[edge[1], partition] <= y[edge])
					model.addConstr( x[edge[1], partition] - x[edge[0], partition] <= y[edge])

					model.addConstr(y[edge] + x[edge[0], partition] + x[edge[1], partition] <= 2 )


			for partition in self.partitions:
				if partition == self.graph.nodes[self.fixed_vertex]["partition"]:
					x[self.fixed_vertex, partition].lb 	= 1
				else:
					x[self.fixed_vertex, partition].ub 	= 0


			if self.Params.Clique_Constraints == True:
				for clique in nx.find_cliques(self.graph):
					clique_list 			= sorted(list(clique))
					clique_size 			= len(clique_list)
					(quotient, reminder)	= divmod(clique_size, self.num_partitions)

					clique_constraint 		= model.addConstr(quicksum(y[edge] for edge in itertools.combinations(clique_list, 2) ) <=  ((clique_size*(clique_size - 1)) - (quotient*(quotient - 1)*(self.num_partitions - reminder) + quotient*(quotient + 1)*reminder) )/2 )
					clique_constraint.Lazy 	= -1

			
			# if self.Params.Wheel_Constraints == True:
			# 	for vertex in self.vertices:
			# 		neighbors 						= list(self.graph.neighbors(vertex))

			# 		induced_graph 					= self.graph.subgraph(neighbors)
			# 		cycle_edges 					= nx.find_cycle(induced_graph)
			# 		largest_cycle 					= cycle_edges[0]
			# 		for edge in cycle_edges:
			# 			largest_cycle += edge

			# 		largest_cycle 					= list(set(largest_cycle))
					
			# 		# bidircted_induced_graph 		= nx.DiGraph(induced_graph)
			# 		# cycles  						= sorted(list(nx.simple_cycles(bidircted_induced_graph)), key = lambda s: len(s))

			# 		# largest_cycle 					= cycles[-1] if cycles else []
			# 		largest_cycle_size 				= len(largest_cycle)

			# 		if largest_cycle_size >= 3:
			# 			cycle_edges 				= zip(largest_cycle, largest_cycle[1:] + [largest_cycle[0]])
			# 			wheel_constraint 			= model.addConstr(quicksum(y[tuple(sorted(edge))] for edge in cycle_edges) - quicksum(y[tuple(sorted((vertex, neighbor) ))] for neighbor in largest_cycle) <=  np.floor(0.5 * largest_cycle_size) )
			# 			wheel_constraint.Lazy 		= -1


			# if self.Params.Bicycle_Wheel_Constraints == True:
			# 	for (vertex1, vertex2) in self.edges:
			# 		common_neighbors 				= list( set(self.graph.neighbors(vertex1)).intersection(set(self.graph.neighbors(vertex2))) )
			# 		induced_graph 					= self.graph.subgraph(common_neighbors)

			# 		bidircted_induced_graph 		= nx.DiGraph(induced_graph)
			# 		cycles  						= sorted(list(nx.simple_cycles(bidircted_induced_graph)), key = lambda s: len(s))

			# 		largest_cycle 					= cycles[-1] if cycles else []
			# 		largest_cycle_size 				= len(largest_cycle)

			# 		if largest_cycle_size >= 3:
			# 			cycle_edges 				= zip(largest_cycle, largest_cycle[1:] + [largest_cycle[0]])
			# 			bicycle_wheel_constraint 	= model.addConstr(quicksum(y[tuple(sorted(edge))] for edge in cycle_edges) + y[(vertex1, vertex2)] - quicksum(y[tuple(sorted((vertex1, neighbor) ))] + y[tuple(sorted((vertex2, neighbor) ))] for neighbor in largest_cycle) <=  2 * np.floor(0.5 * largest_cycle_size) - largest_cycle_size + 1)
			# 			bicycle_wheel_constraint.Lazy 	= -1



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

			model 							= model if self.Params.Relaxed == False else model.relax()
			model.update()


			#-----------------------------------------------------------------------------------
			# Solve the model and extract the obtained solution
			#-----------------------------------------------------------------------------------
			model.optimize()

			self.gurobi_obj_value			= model.objVal
			self.gurobi_BB_nodes 			= model.getAttr(GRB.Attr.NodeCount) 		# Spatial BB nodes
			self.gurobi_ObjBound			= model.ObjBound
			self.gurobi_MIPGap				= model.MIPGap if self.Params.Relaxed == False else 0.0


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
					
	
	
	






