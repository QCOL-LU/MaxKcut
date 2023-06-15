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

		self.Twin_Fix						= False
		self.Edge_Based_Fix					= False
		self.Rehfeldt_Fix 					= False

		self.Lange_Fold 					= False


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
		self.QAOA_Num_Shots 				= 10000
		self.QAOA_Scipy_Optimizer			= "COBYLA"
		self.QAOA_Opt_Print_Time			= 1
		self.QAOA_Opt_Tol 					= 1e-3
		self.QAOA_Verbosity 				= 1
		self.QAOA_Angles					= []
		self.QAOA_Brute_Num_Samples 		= 50
		self.QAOA_Gates_Noise 				= 0.0 	# Error probabilities; 0.01 means that the error is 1%
		self.Penalty_Increase 				= 1
		self.Naive 							= False
		self.Adjusted_Penalty				= False

		self.Grad_Epoch_Len 				= 20
		self.Grad_Method 					= "RandGradEst"
		self.Grad_Sample_Size 				= 20
		self.Grad_Smoothing_Param 			= 1e-1



		self.Gurobi_TimeLimit 				= 3600
		self.Gurobi_MIPGap 					= 0
		self.Gurobi_Threads					= 10
		self.Gurobi_Heuristics 				= 0.05
		self.Gurobi_Presolve 				= -1
		self.Gurobi_Symmetry				= -1
		self.Gurobi_Cuts 					= -1
		self.Gurobi_LogToConsole 			= 0


		self.Is_Simulator 					= True
		self.QLSA_First_Run 				= True
		self.Hub 							= 'ibm-q-ncsu' 			# Confidential *
		self.Group 							= 'lehigh-universit'	# Confidential *
		self.Project 						= 'qc-for-comb-opt'		# Confidential *
		self.QC_Name 						= 'ibm_washington' #'ibmq_qasm_simulator'		#'ibmq_montreal'	ibm_washington	ibm_sherbrooke # Confidential * 
		self.Token 							= "cb03fdfa5944bde5276681dba3d7c581b603da237f410af03524ea27fc94cb97b1d98c4190b66808382dcf517015bc36968dc0f4426783889390c8265701f6cc"

