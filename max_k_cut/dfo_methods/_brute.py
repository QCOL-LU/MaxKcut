import time
import itertools
from copy import deepcopy

import numpy as np


#================================================================================================
# The Partition-based MILO formulation
#================================================================================================
def optimize_brute_search(self):

	bounds 					= [(0,  2*np.pi) for _ in range(self.Params.QAOA_Num_Levels)] + [(0,  np.pi) for _ in range(self.Params.QAOA_Num_Levels)]
	Ns						= self.Params.QAOA_Brute_Num_Samples



	# optimizer_results 		= brute(func=self.qaoa_expected_value, 
	# 								ranges=bounds,
	# 								Ns=10)

	self.dfo_start_time 	= time.time()

	#--------------------------------------------------------------------------------------------
	# Model initialization
	#--------------------------------------------------------------------------------------------
	iter_finish				= 100

	if self.Params.QAOA_Num_Levels == 1:
		betas, gammas		= np.meshgrid(np.linspace(bounds[1][0]/np.pi, bounds[1][1]/np.pi, int(Ns/2)), np.linspace(bounds[0 ][0]/np.pi, bounds[0][1]/np.pi, Ns ))
		
		objectives 			= np.zeros((Ns, int(Ns/2)))


	itr 					= 0
	total_itr 				= 0
	for g_ind, gamma1 in enumerate(np.linspace(*bounds[0], Ns) ):
		for b_ind, beta1 in enumerate(np.linspace(*bounds[self.Params.QAOA_Num_Levels ], int(Ns/2))):

			if self.Params.QAOA_Num_Levels == 1:
				old_best_exp_obj			= self.qaoa_best_avg_obj_value

				arg 						= (gamma1, beta1)
				obj 						= - self.qaoa_expected_value(arg)
				total_itr 					+= 1
				objectives[g_ind][b_ind]	= obj

				if self.qaoa_best_avg_obj_value != old_best_exp_obj:
					itr 					= 0
				else:
					itr 					+= 1


			if self.Params.QAOA_Num_Levels >= 2:
				for gamma2 in np.linspace(*bounds[1], Ns):
					for beta2 in np.linspace(*bounds[self.Params.QAOA_Num_Levels + 1], int(Ns/2)):

						if self.Params.QAOA_Num_Levels == 2:
							arg 	= (gamma1, gamma2, beta1, beta2)
							obj 	= - self.qaoa_expected_value(arg)

						if self.Params.QAOA_Num_Levels == 3:
							for gamma3 in np.linspace(*bounds[2], Ns):
								for beta3 in np.linspace(*bounds[self.Params.QAOA_Num_Levels + 2], int(Ns/2)):

								
									arg 	= (gamma1, gamma2, gamma3, beta1, beta2, beta3)
									obj 	= - self.qaoa_expected_value(arg)

			
	results  						= type('', (), {})()

	results.gammas 					= gammas
	results.betas 					= betas

	results.objectives				= objectives
	results.locally_optimal_angles 	=  self.best_angles
	results.message 				= "The brute search stopped after "+ str(total_itr) + " iterations."

	self.plot_qaoa_level_one(gammas, betas, objectives)

	self.dfo_end_time 	= time.time()

	return results
					
	

