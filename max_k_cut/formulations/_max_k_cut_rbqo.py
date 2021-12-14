import networkx as nx 
from gurobipy import *
import time
import itertools


#================================================================================================
# Solve the R-BQO formulation
#================================================================================================
def solve_max_k_cut_rbqo(self):

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
			x 		= model.addVars(self.vertices, self.partitions[:-1], vtype=GRB.BINARY, name="x")

			
			#-----------------------------------------------------------------------------------
			# Objective function
			#-----------------------------------------------------------------------------------
				
			if self.Params.Curvature_Type == "indefinite":

				model.setObjective(quicksum( self.graph.edges[edge]["weight"] * (quicksum( x[edge[0], partition] + x[edge[1], partition] - 2*x[edge[0], partition] * x[edge[1], partition] for partition in self.partitions[:-1]) \
									- quicksum(x[edge[0], partition1] * x[edge[1], partition2] + x[edge[0], partition2] * x[edge[1], partition1] for (partition1, partition2) in itertools.combinations(self.partitions[:-1], 2)) )
									for edge in self.edges), GRB.MAXIMIZE)		
			
			elif self.Params.Curvature_Type in ["convex", "concave"]:

				self.calculate_curvature_coefs()

				model.setObjective(quicksum( self.graph.edges[edge]["weight"] * (quicksum( x[edge[0], partition] + x[edge[1], partition] - 2*x[edge[0], partition] * x[edge[1], partition] for partition in self.partitions[:-1]) \
											- quicksum(x[edge[0], partition1] * x[edge[1], partition2] + x[edge[0], partition2] * x[edge[1], partition1] for (partition1, partition2) in itertools.combinations(self.partitions[:-1], 2)) ) for edge in self.edges) \
											+ quicksum(self.variable_coefs[vertex][partition] * (x[vertex, partition] *  x[vertex, partition] - x[vertex, partition]) for vertex in self.vertices for partition in self.partitions[:-1]), \
											GRB.MAXIMIZE)
			else:
				sys.exit("Please specify the correct objective type: (indefinite, convex)")
			
			
			model.update()
			#-----------------------------------------------------------------------------------
			# Constraints
			#-----------------------------------------------------------------------------------
			for vertex in self.vertices:
				model.addConstr(quicksum(x[vertex, partition] for partition in self.partitions[:-1]) <= 1)

				if self.Params.Relaxed == True and self.Params.Relaxed_NonConvex == True:
					for partition in self.partitions[:-1]:
						model.addConstr(x[vertex][partition] <= x[vertex][partition]*x[vertex][partition])


			model.update()
			#-----------------------------------------------------------------------------------
			# Set the solver parameters
			#-----------------------------------------------------------------------------------
			model.Params.Heuristics 		= self.Params.Gurobi_Heuristics
			model.Params.Presolve 			= self.Params.Gurobi_Presolve
			model.Params.Symmetry 			= self.Params.Gurobi_Symmetry
			model.Params.Cuts 				= self.Params.Gurobi_Cuts

			model.Params.NonConvex 			= -1 if self.Params.Curvature_Type == "convex" else 2
			
			model.Params.Threads			= self.Params.Gurobi_Threads
			model.Params.timeLimit 			= self.Params.Gurobi_TimeLimit
			model.Params.MIPGap 			= self.Params.Gurobi_MIPGap

			# model.Params.LogToConsole 		= 0
			model.Params.LogFile 			= self.filename[:-4] + "_log.txt"
			
			model 							= model.relax() if self.Params.Relaxed == True else model
			model.update()


			#-----------------------------------------------------------------------------------
			# Solve the model and extract the obtained solution
			#-----------------------------------------------------------------------------------
			model.optimize()

			self.gurobi_obj_value			= model.objVal
			self.gurobi_BB_nodes 			= model.getAttr(GRB.Attr.NodeCount) 		# Spatial B&B nodes
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
			solution						= {vertex: [model.getVarByName("x["+str(vertex)+","+str(partition)+"]").x for partition in self.partitions[:-1]] for vertex in self.vertices}
	
	
	for vertex in self.vertices:
		solution[vertex].append(1 - sum(solution[vertex]) )

	if self.Params.Rounding_Heuristic == True and self.Params.Relaxed == True: 
		for vertex in self.vertices:
			b_vector 						= [sum([self.graph.edges[(vertex, neighbor)]["weight"] * solution[neighbor][partition] for neighbor in self.graph.neighbors(vertex)]) for partition in self.partitions]
			selected_partition 				= b_vector.index(min(b_vector))
			solution[vertex] 				= [1 if partition == selected_partition else 0 for partition in self.partitions]
			

	if self.Params.Rounding_Heuristic == True or self.Params.Relaxed == False:

		for vertex in self.vertices:
			for partition in self.partitions:

				if solution[vertex][partition] > 0.5:
					self.graph.nodes[vertex]["partition"] 	= partition 


	self.gurobi_end_time	 		= time.time()

	self.print_gurobi_results_summary()
					
	


			

	