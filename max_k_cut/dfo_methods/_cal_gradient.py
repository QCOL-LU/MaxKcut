import time
import itertools
from copy import deepcopy

import numpy as np

from max_k_cut.quantum_methods._max_k_cut_qaoa_circuits import qaoa_expected_value

#================================================================================================
# Calculate the estimated gradient
#================================================================================================
def cal_gradient(self, angles):
	dimension 				= 2 * self.Params.QAOA_Num_Levels
	smoothing_param 		= self.Params.Grad_Smoothing_Param
	
	if self.Params.Grad_Method == "RandGradEst":
		np.random.seed(seed=self.seed)
		random_direction	= np.random.rand(dimension)
		random_direction 	= random_direction / np.linalg.norm(random_direction)
		
		f_angle 			= self.qaoa_expected_value(angles)
		f_new_angle			= self.qaoa_expected_value(angles + smoothing_param * random_direction)

		gradient 			= (dimension / smoothing_param) * (f_new_angle - f_angle) * random_direction

	elif self.Params.Grad_Method == "Avg-RandGradEst":
		sample_size 		= self.Params.Grad_Sample_Size
		gradient_sum 		= np.zeros(dimension)

		for _ in range(sample_size):
			self.Params.Grad_Method 	= "RandGradEst"
			gradient_sum 	+= self.cal_gradient(angles)
			self.seed 		+= 1

		self.Params.Grad_Method 		= "Avg-RandGradEst"
		gradient 			= (dimension / (smoothing_param * sample_size) ) * gradient_sum


	elif self.Params.Grad_Method == "CoordGradEst":

		gradient_sum 		= np.zeros(dimension)

		for i in range(dimension):
			direction 		= np.zeros(dimension)
			direction[i] 	= 1

			f_angle_minus 	= self.qaoa_expected_value(angles + smoothing_param * direction)
			f_angle_plus 	= self.qaoa_expected_value(angles - smoothing_param * direction)
			gradient_sum	+= (f_angle_plus - f_angle_minus) * direction

		gradient	 		= (1 / (2 * smoothing_param)) * gradient_sum



	return gradient
					
	
#================================================================================================
# Calculate the estimated gradient of a sample batch
#================================================================================================
def cal_gradient_batch(self, angles, batch_size):
	dimension 				= 2 * self.Params.QAOA_Num_Levels
	gradient_sum 			= np.zeros(dimension)

	for i in range(batch_size):
		gradient_sum 		+= self.cal_gradient(angles)


	gradient 				= (1 / batch_size) * gradient_sum

	return gradient 






