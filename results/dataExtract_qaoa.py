import pandas as pd
import os.path
from os import walk
import os 
from copy import deepcopy

files 				= []
wrong_words 		= ["log", ".py"]

for folder in sorted(os.listdir(os.getcwd())):
	if not os.path.isdir(folder) or "k" not in folder: continue
	folder_paths 	= os.path.join(os.getcwd(), folder)	


	for (ind, filename) in enumerate(sorted(os.listdir(folder_paths))):

		# if "preprocess" in filename:
		# 	os.remove(folder_paths+"/"+filename)
		# 	continue
		
		flag 		= 0
		for word in wrong_words:
			if word in filename: flag = 1
		if flag == 1: continue

		file_address = os.path.join(folder_paths, filename)

		if ((not os.path.isdir(file_address)) and not ("log" in filename)):
			files.append(file_address )


# headers12	= ["File name", "Name of instance", "Density (%)","Naive Penalty Coef", "Is planar", "Number of vertices", "Number of edges", "Number of partitions", "Method","Peel", "Decompose", "Fold", "Penalty Multiplier", "Gates Error Probability",
# 				"Avg of QAOA obj (BQO obj)", "Avg of QAOA obj (penalty)", "Avg of constraint violation","Feasible percentage (%)", "Avg of pure feasible obj", 
# 				"Avg of QAOA obj", "Avg of QAOA obj (feasible)", "STD of QAOA obj", "STD of QAOA obj (feasible)", "Best QAOA obj", "Modified best QAOA obj", "Best QAOA obj (feasible)", "QAOA total time",  "QAOA circuit depth", "Num of bi-comp", "Vertex num in largest comp", "Edge num in largest comp",  "Pre-processing running time", "Total solver time", "Running time", "Upper bound" ,"Objective value"]

headers12	= ["File name", "Name of instance", "Density (%)","Naive Penalty Coef",  "Number of vertices", "Number of edges", "Number of partitions", "Method",
				"Avg of QAOA obj (BQO obj)", "Avg of QAOA obj (penalty)", "Avg of constraint violation","Feasible percentage (%)", "Avg of pure feasible obj", 
				"Avg of QAOA obj", "Avg of QAOA obj (feasible)", "STD of QAOA obj", "STD of QAOA obj (feasible)", "Best QAOA obj", "Modified best QAOA obj", 
				"Avg of tight R-QUBO QAOA obj (penalty)", "Avg of tight R-QUBO QAOA obj", "Avg of tight QAOA obj (penalty)", "Avg of tight QAOA obj",
				"Best QAOA obj (feasible)", "QAOA total time",  "QAOA circuit depth",  "Total solver time", "Running time", "Upper bound" ,"Objective value"]



string_title 		= deepcopy(headers12) + ["Method","Peel", "Decompose", "Fold", "QAOA Num Levels", "QAOA Scipy Optimizer", "Gamma angles", "Beta angles", "Naive Penalty Coef"]
large_headers 		= deepcopy(headers12) 


empty_rows 			= [["" for item in large_headers] for row in range(10000)]
df 					= pd.DataFrame(empty_rows, columns=large_headers)




for (row, filename) in enumerate(files):

	with open(filename) as f:
		file_name 			= filename.split("/")[-1]
		df.iloc[row , df.columns.get_loc("File name")] 	= file_name
		# print(file_name)
		lines 				= f.read().splitlines()
		lines.reverse()
		
		for header in headers12:
			for line in lines:
				title 		= line.split(":")[0]


	
				if title == header:
					
					if title in string_title:
						value 				= line.split(":")[1].strip()

					else:
						value 				= float(line.split(":")[1].strip())

					df.iloc[row, df.columns.get_loc(title)] 	= value 
					break


df.to_csv('results.csv')			


