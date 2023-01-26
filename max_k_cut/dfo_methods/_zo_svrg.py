import time
import itertools
from copy import deepcopy

import numpy as np
from max_k_cut.dfo_methods._cal_gradient import *
from max_k_cut.dfo_methods._dfo_results import *
import random



#================================================================================================
# The Partition-based MILO formulation
#================================================================================================
def back_to_range(angles):
	new_angles 		= [angle - np.floor(angle/(2 * np.pi)) * 2 * np.pi for angle in angles]
	return new_angles

#================================================================================================
# The Partition-based MILO formulation
#================================================================================================
def optimize_zo_svrg(self):

	optimizer_results  				= Optimizer()

	self.dfo_start_time 			= time.time()

	angles 							= np.ones(2 * self.Params.QAOA_Num_Levels) * np.pi / 2
	step_size 						= 5e-2 

	init_QAOA_Num_Shots				= self.Params.QAOA_Num_Shots
			
	total_size 						= 10
	self.Params.QAOA_Num_Shots 		= 500


	batch_size 						= init_QAOA_Num_Shots/(self.Params.QAOA_Num_Shots * total_size)
	gradient 						= self.cal_gradient_batch(angles, total_size)

	next_iterates					= [angles for _ in range(self.Params.Grad_Epoch_Len)]


	while  (step_size <= 1e-5): #(p.linalg.norm(gradient) > self.Params.QAOA_Opt_Tol):
		
		gradient_norm_init 			= np.linalg.norm(gradient)

		gradient_batch_init 		= self.cal_gradient_batch(angles, batch_size)
		angles_epoch				= deepcopy(angles)

		for i in range(self.Params.Grad_Epoch_Len):
			zo_gradient_blending 	= self.cal_gradient_batch(angles, batch_size) - gradient_batch_init + gradient
			angles_epoch 			= back_to_range(angles_epoch + step_size * zo_gradient_blending)


			gradient_norm 				= np.linalg.norm(zo_gradient_blending)
			gradient_ratio 				= gradient_norm/gradient_norm_init

			step_size 					= max((max(.5, gradient_ratio) if gradient_ratio < .3 else 1.5) * step_size, 1e-5)




		angles 						= deepcopy(back_to_range(angles_epoch) )
		gradient 					= self.cal_gradient_batch(angles, total_size)


		

		

	rand_index 						= random.randint(0, self.Params.Grad_Epoch_Len)

	optimizer_results.angles 		= next_iterates[rand_index]


	self.Params.QAOA_Num_Shots 		= init_QAOA_Num_Shots
	
	self.dfo_end_time 				= time.time()

	optimizer_results.run_time 		= self.dfo_end_time - self.dfo_start_time

	optimizer_results.message 		= "Optimizer stopped after "+ str(self.qaoa_iter) + " iterations.\n"

	return optimizer_results
					
	

