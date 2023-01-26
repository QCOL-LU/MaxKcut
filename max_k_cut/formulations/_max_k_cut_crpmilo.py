import networkx as nx 
from gurobipy import *
import time
import itertools
from copy import deepcopy
import random


from networkx.algorithms import tree

#================================================================================================
# The Reduced Partition-based MILO formulation proposed by Wang and Hijazi
#================================================================================================
def solve_max_k_cut_crpmilo(self):

	self.gurobi_start_time 			= time.time()

	#--------------------------------------------------------------------------------------------
	# Construct the extended chordal graph
	#--------------------------------------------------------------------------------------------
	chordal_graph 					= deepcopy(self.graph)
	chordal_graph_edges 			= deepcopy(self.edges)

	if not nx.is_chordal(chordal_graph):
		temp_graph					= deepcopy(self.graph)

		while (temp_graph.number_of_nodes() >= 2):
			vertex_neighbor			= {vertex: list(temp_graph.neighbors(vertex)) for vertex in temp_graph.nodes}
			num_fill_in_edges 		= {vertex: len(vertex_neighbor[vertex]) * (len(vertex_neighbor[vertex]) - 1)/2 - len(list(temp_graph.subgraph(vertex_neighbor[vertex]).edges)) for vertex in temp_graph.nodes}

			vertex  				= min(num_fill_in_edges, key=num_fill_in_edges.get)

			
			if num_fill_in_edges[vertex] > 0:

				for (neighbor1, neighbor2) in itertools.combinations(vertex_neighbor[vertex], 2):
					chordal_graph.add_edge(min(neighbor1, neighbor2), max(neighbor1, neighbor2) )
					temp_graph.add_edge(min(neighbor1, neighbor2), max(neighbor1, neighbor2) )
					chordal_graph_edges.append((min(neighbor1, neighbor2), max(neighbor1, neighbor2)))

			
			temp_graph.remove_node(vertex)

			if max(num_fill_in_edges.values()) == 0: break


	maximal_cliques 					= list(nx.find_cliques(chordal_graph))
	
	self.density_chordal_graph 			= 200 * chordal_graph.number_of_edges() / ( (self.num_vertices - 1) * self.num_vertices )
	directed_graph 						= self.graph.to_directed()


	spanning_tree_edge		= list(tree.minimum_spanning_edges(chordal_graph, algorithm="kruskal", data=False) )
	forest_edges			= random.sample(spanning_tree_edge, self.num_vertices - self.num_partitions)
	forest 					= nx.Graph()
	forest.add_edges_from(forest_edges)

	components 				= list(nx.connected_components(forest))



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
							for (vertex1, vertex2) in chordal_graph_edges}

			t 			= {edge: model.addVar(vtype=GRB.BINARY, name="z(%i,%i)" %edge) \
							for edge in directed_graph.edges()}


			temp_edges 	= []

			for component in components:
				for (vertex1, vertex2) in itertools.combinations(component, 2):
					if (vertex1, vertex2) in chordal_graph_edges:
						z[(vertex1, vertex2)].start = 1
						temp_edges.append((vertex1, vertex2))
						# print((vertex1, vertex2))


					if (vertex2, vertex1) in chordal_graph_edges:
						z[(vertex2, vertex1)].start = 1 
						temp_edges.append((vertex2, vertex1))
						# print((vertex1, vertex2))


			for edge in chordal_graph_edges:
				if edge not in temp_edges:
					z[edge].start = 0
					if edge in self.edges:
						t[(edge[0], edge[1])].start = 0
						t[(edge[1], edge[0])].start = 0




			model._z 	= z
			model._t 	= t

			#-----------------------------------------------------------------------------------
			# Objective function
			#-----------------------------------------------------------------------------------
			model.setObjective(quicksum(self.graph.edges[edge]["weight"] * (1 - z[edge]) for edge in self.edges), GRB.MAXIMIZE)
			

			#-----------------------------------------------------------------------------------
			# Constraints
			#-----------------------------------------------------------------------------------
			for clique in maximal_cliques:
				if len(clique) >= 3:
					for vertex_set in itertools.combinations(clique, 3):
						vertex_set 			= sorted(vertex_set)
						model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[0], vertex_set[2])] <= 1 + z[(vertex_set[1], vertex_set[2])])
						model.addConstr(z[(vertex_set[0], vertex_set[1])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[2])])
						model.addConstr(z[(vertex_set[0], vertex_set[2])] + z[(vertex_set[1], vertex_set[2])] <= 1 + z[(vertex_set[0], vertex_set[1])])


			model.addConstrs(t[(vertex1, vertex2)] + t[(vertex2, vertex1)] <=  z[(vertex1, vertex2)] for (vertex1, vertex2) in self.edges)

			model.addConstr(quicksum(t[edge] for edge in directed_graph.edges()) == self.num_vertices - self.num_partitions )

			model.addConstrs(quicksum(t[(neighbor, vertex)] for neighbor in self.graph.neighbors(vertex)) <= 1 for vertex in self.vertices)





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

					# for clique in maximal_cliques:
					# 	if len(clique) <= self.num_partitions:
					# 		continue
							
					# 	clique.sort()
					# 	z_star_edges 		= [(vertex1, vertex2) for (vertex1, vertex2) in itertools.combinations(clique, 2) \
					# 							 if (model.cbGetSolution(model._z[(vertex1, vertex2)]) ) > .5]

					# 	graph 				= nx.Graph()

					# 	graph.add_nodes_from(clique)
					# 	graph.add_edges_from(z_star_edges)

					# 	components 			= [list(component)[0] for component in nx.connected_components(graph)]

					# 	size 				= len(components)

					# 	if size <= self.num_partitions:
					# 		continue


					# 	for (ind, _) in enumerate(components):
					# 		first 			= components[ind: min(size, ind + self.num_partitions + 1)]
					# 		size_second 	= self.num_partitions + 1 - len(first)
					# 		second 			= components[:size_second]
					# 		set_Q 			= first + second

					# 		model.cbLazy(quicksum(model._z[(min(key1, key2), max(key1, key2) )] \
					# 				for (ind1, key1) in enumerate(set_Q[:-1]) for key2 in set_Q[ind1 + 1:] ) >= 1)
					# 		model.update()



					t_star_edges 		= [(vertex1, vertex2) for (vertex1, vertex2) in directed_graph.edges() \
													 if (model.cbGetSolution(model._t[(vertex1, vertex2)]) ) > .5]

					graph_subtour 		= nx.Graph()
					graph_subtour.add_nodes_from(self.vertices)
					graph_subtour.add_edges_from(t_star_edges)

					components 			= list(nx.connected_components(graph_subtour))

					if len(components) != self.num_partitions: 

						for component in components:
							component_edges = list(graph_subtour.subgraph(component).edges())

							if len(component_edges) >= len(component):
								model.cbLazy(quicksum(model._t[(vertex1, vertex2)] + model._t[(vertex2, vertex1)] \
											for (vertex1, vertex2) in component_edges) <= len(component) - 1)
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

			model.Params.symmetry = 2
			model.update()

			#-----------------------------------------------------------------------------------
			# Solve the model and extract the obtained solution
			#-----------------------------------------------------------------------------------
			# model.optimize(add_lazy_constraint)
			# model.optimize()
			model.optimize(add_lazy_constraint)

			self.gurobi_obj_value			= model.objVal
			self.gurobi_BB_nodes 			= model.getAttr(GRB.Attr.NodeCount) 		# Spatial BB nodes
			self.gurobi_ObjBound			= model.ObjBound
			self.gurobi_MIPGap				= model.MIPGap
			
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
			z_star_edges 					= [(vertex1, vertex2) for (vertex1, vertex2) in chordal_graph_edges if (z[(vertex1, vertex2)].x) > 0.5]

			conflict_vertex_pairs 			= [(vertex1, vertex2) for (vertex1, vertex2) in chordal_graph_edges if (z[(vertex1, vertex2)].x) <= 0.5]

			graph 							= nx.Graph()

			graph.add_nodes_from(self.vertices)
			graph.add_edges_from(z_star_edges)

			components 						= [list(component) for component in nx.connected_components(graph)]
			num_components 					= len(components)

			if num_components == self.num_partitions:
				print("connected")
			else:
				print("not connected", num_components)

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
					
	

