import xlwt 
import os.path
from os import walk
import os 

files 				= []
wrong_words 		= ["log", ".py"]

for folder in sorted(os.listdir(os.getcwd())):
	if not os.path.isdir(folder) or "max" not in folder: continue
	folder_paths 	= os.path.join(os.getcwd(), folder)	

	for (ind, filename) in enumerate(sorted(os.listdir(folder_paths))):
		flag 		= 0
		for word in wrong_words:
			if word in filename: flag = 1
		if flag == 1: continue

		if (not os.path.isdir(filename) and not ("log" in filename)):
			files.append(os.path.join(folder_paths, filename) )

book 		= xlwt.Workbook(encoding="utf-8")
sheet1 		= book.add_sheet("results", cell_overwrite_ok=True)
sheet2 		= book.add_sheet("models", cell_overwrite_ok=True)
headers12	= {"Instance parameters": ["Name of instance", "Number of nodes", "Number of edges", "Number of partitions", "Category"], "Solver parameters":["Method","Peel", "Decompose", "Fold", "Curvature Type", "Curvature Method", "Clique Constraints"]}

headers1 	= ["Name", "Vertex num in largest comp", "Edge num in largest comp",  "Pre-processing Running time", "Running time", "Objective value"] 
headers2 	= ["Name", "Model status", "Explored B&B nodes", "Objective value", "Upper bound", "Optimality gap", "Gurobi running time", ]


sheet1.write_merge(0,1, 0,0, "File name")
sheet2.write_merge(0,1, 0,0, "File name")
sheet1.write_merge(0,0, 1,4, "Instance parameters")
sheet2.write_merge(0,0, 1,4, "Instance parameters")
sheet1.write_merge(0,0, 5,11, "Solver parameters")
sheet2.write_merge(0,0, 5,11, "Solver parameters")

for (ind, header) in enumerate(headers12["Instance parameters"]):
	sheet1.write(1, ind + 1, header)
	sheet2.write(1, ind + 1, header)

for (ind, header) in enumerate(headers12["Solver parameters"]):
	sheet1.write(1, ind + 5, header)
	sheet2.write(1, ind + 5, header)


def is_header(line, cols, col_heads, itr, main_header):
	header_changed			= [False, False]

	if col_heads[0] == 0 or col_heads[1] == 0:
		sheet_list 			= []
		col_heads			= [1, 1]

	elif (itr > 0 and "====" in line) or ("0" not in line):
		itr 				= -1
		sheet_list 			= []
		main_header			= ""
	
	elif "Instance parameters" in line:
		header_changed		= [True, True]
		sheet_list 			= [0, 1]
		itr 				= 4
		main_header			= "Instance parameters"
		col_heads			= [1, 1]


	
	elif "Solver parameters" in line:
		header_changed		= [True, True]
		sheet_list 			= [0, 1]
		itr 				= 7
		main_header			= "Solver parameters"
		col_heads			= [11, 11]

	
	elif "Summary of results" in line:

		if "model" in line:
			header_changed[1] 	= True
			sheet_list 		= [1]
			main_header		= "Summary of results of model"
			itr 			= 7

			sheet2.write_merge(0,0, col_heads[1]+1,col_heads[1]+7, "Summary of results of model")

			for (ind, header) in enumerate(headers2):
				sheet2.write(1, ind + col_heads[1] + 1, header)

			col_heads[1] 		+= 7
				


		elif "QAOA" in line: 
			print("add QAOA")

		
		else:
			header_changed[0] 	= True
			sheet_list 			= [0]
			main_header			= "Summary of results"
			itr 				= 5
			sheet1.write_merge(0,0, col_heads[0]+1,col_heads[0]+5, "Summary of results")

			for (ind, header) in enumerate(headers1):
				sheet1.write(1, ind + col_heads[0] + 1, header)

			col_heads[0] 		+= 5
	else:
		itr 				= -1
		sheet_list 			= []
		main_header			= ""
	return sheet_list, itr, col_heads, main_header, header_changed





for (row, filename) in enumerate(files):
	with open(filename) as f:

		file_name 	= filename.split("/")[-1]
		sheet1.write(row+2, 0, file_name)
		sheet2.write(row+2, 0, file_name)
		main_header	= ""
		
		lines 		= f.read().splitlines()

		itr 		= -1
		
		cols 		= [0, 0]

		col_heads	= [0, 0]


		for line in lines:

			(sheet_list, itr, col_heads, main_header, header_changed)		= is_header(line, cols, col_heads, itr, main_header)
			print(row, col_heads, filename)
			
			title 		= line.split(":")[0]


			if itr > -1 :
				
				if main_header == "Instance parameters":
					header_list 	= headers12[main_header]
					ind 			= 1

				elif main_header == "Solver parameters":
					header_list 	= headers12[main_header]
					ind 			=  5

				elif main_header == "Summary of results":
					header_list 	= headers1
					ind 			= col_heads[0] - 5


				elif main_header == "Summary of results of model":
					header_list 	= headers2
					ind 			= col_heads[1] - 7

				if title in header_list or "Summary" in title:
					col 		= ind

					if title in ["Name of instance", "Model status"] or main_header  == "Solver parameters":
						val 	= line.split(":")[0].strip()
						col 	+= header_list.index(title) 
					
					elif "Summary" in title and "model" in title:
						val 	= line.split("of")[-1].strip()

					elif "Summary" in title:
						val 	= line.split("of")[-1].strip()
					
					else:
						val 	= float(line.split(":")[0].strip())
						col 	+= header_list.index(title) 


					if 0 in sheet_list:
						sheet1.write(row+2, col, val)
						cols[0] 		= col

					if 1 in sheet_list:
						sheet2.write(row+2, col, val)
						cols[1] 		= col


					itr 		-= 1
			
book.save("results.xls")	


