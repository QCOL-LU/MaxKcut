from copy import deepcopy

import numpy as np
import random


#===========================================================================
# Linear optimization problem class 
#===========================================================================
class Optimizer():

	#-----------------------------------------------------------------------
	# Initialize the Optimizer results
	#-----------------------------------------------------------------------
	def __init__(self):

		self.message 	= ""
		self.angles 	= []
		self.run_time 	= 0.0