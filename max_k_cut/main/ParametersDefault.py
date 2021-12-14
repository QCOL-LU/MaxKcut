#===========================================================================
# The class of parameters 
#===========================================================================
class Parameters():

	#=======================================================================
	# Class initialization
	#=======================================================================
	def __init__(self):

		#-------------------------------------------------------------------
		# Default values of paramters
		#-------------------------------------------------------------------
		self.Method 						= "BQO"
		self.Peel 							= True
		self.Decompose 						= True
		self.Fold 							= True
		self.Bipartite 						= True


		self.Num_Partitions 				= 2
		self.Verbosity 						= 1


		self.Relaxed 						= False
		self.Relaxed_NonConvex 				= False
		self.Rounding_Heuristic 			= False

		self.Curvature_Type 				= "indefinite"		# indefinite 	convex			concave 
		self.Curvature_Method 				= "power_iteration"	# mosek		power_iteration		numpy
		
		self.Symmetry_Breaking 				= False
		self.Clique_Constraints 			= True
		self.Wheel_Constraints 				= True
		self.Bicycle_Wheel_Constraints		= True

		self.QAOA_Optimize 					= True
		self.QAOA_Num_Levels 				= 3
		self.QAOA_Num_Shots 				= 1000
		self.QAOA_Scipy_Optimizer			= "COBYLA"
		self.QAOA_Opt_Print_Time			= 1
		self.QAOA_Opt_Tol 					= 1e-3
		self.QAOA_Verbosity 				= 1
		self.QAOA_Angles					= []
		self.Penalty_Increase 				= 1


		self.Gurobi_TimeLimit 				= 3600
		self.Gurobi_MIPGap 					= 0
		self.Gurobi_Threads					= 10
		self.Gurobi_Heuristics 				= 0.05
		self.Gurobi_Presolve 				= -1
		self.Gurobi_Symmetry				= -1
		self.Gurobi_Cuts 					= -1
		self.Gurobi_LogToConsole 			= 0





