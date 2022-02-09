import networkx as nx 
import numpy as np
import sys 


from copy import deepcopy

from qiskit.visualization import circuit_drawer
import matplotlib.pyplot as plt



#================================================================================================
# Plot the graph of objective value vs probability
#================================================================================================
def plot_qaoa_solutions_dist(self, feasible_sol=False, num_bins=20):
	histogram 			= self.best_histogram if feasible_sol==False else self.feasible_best_histogram

	(lb_bin, ub_bin)	= (min(histogram.keys()), max(histogram.keys()))

	range_bin 			= (ub_bin - lb_bin)
	width_bin			= range_bin * (num_bins + 1) / num_bins**2
	precision_bin		= max(-int(np.floor(np.log(width_bin))), 2) 

	(lb_bin, ub_bin)	= (lb_bin - width_bin/2, ub_bin + width_bin/2)
	bounds_bins 		= [ ((lb_bin + width_bin * itr), (lb_bin + width_bin * (itr + 1))) for itr in range(num_bins) ]

	dict_bin			= {}

	for (key, value) in histogram.items():
		for (lb, ub) in bounds_bins:
			if lb <= key and key < ub and value > 0:
				mid_bin				= np.round((lb + ub)/2, precision_bin) 
				dict_bin[mid_bin] 	= dict_bin.get(mid_bin, 0) + value


	min_value 			= max(dict_bin.values()) * .001

	for (key, value) in dict_bin.items():
		if value < min_value:
			dict_bin[key] 			= min_value
	
	plt.bar(dict_bin.keys(), 
			dict_bin.values(), 
			align="center", 
			width=0.85*min([ub - lb for (lb, ub) in bounds_bins]), 
			edgecolor="black", 
			color="lightblue")
	

	plt.ylabel('Probability (%)')
	temp 				= self.qaoa_best_avg_obj_value if feasible_sol==False else self.qaoa_feasible_best_avg_obj_value 
	plt.xlabel('Objective value with EV = ' + str(temp))

	temp 				= '_obj_prob.png' if feasible_sol==False else '_obj_feasible_prob.png'
	plt.savefig(self.figure_path + "/" + self.name + "_" + self.directory+ "_" +self.Params.Method + temp, dpi=300)
	plt.clf()

	filename 			= self.figure_path + "/" + self.name + "_"  + self.directory+ "_" +self.Params.Method +'_circuit'
	# circuit_drawer(self.qaoa_circuit, filename=filename, output='mpl', style={'backgroundcolor': '#EEEEEE'})





#================================================================================================
# Plot the graph 
#================================================================================================
def plot_graph_problem(self):

	default_axes 	= plt.axes(frameon=True)
	pos				= nx.circular_layout(self.graph)

	options 		= {"font_size": 14,
						"node_size": 1000,
						"node_color": "white",
						"edgecolors": "black",
						"linewidths": 1.5,
						"width": 1.5,
						}

	nx.draw_networkx(self.graph, ax=default_axes, pos=pos, **options)

	plt.savefig(self.figure_path + "/"+ self.name + ".png", dpi=300)
	plt.clf()



#================================================================================================
# Plot the figure of the graph corresponding to the obtained solution
#================================================================================================
def plot_graph_solution(self):

	colors 			= ["#CCFFFF", "#E5FFCC", "#CD1076","#FFCCE5", "#FFFF99", "#FFF8DC","#FFE5CC", "#FFCCCC", "#C1CDCD","#FFFF99"]

	node_colors		= [colors[self.graph.nodes[vertex]["partition"]] for vertex in self.vertices]

	not_cut_edges 	= [edge for edge in self.edges\
						if self.graph.nodes[edge[0]]["partition"] == self.graph.nodes[edge[1]]["partition"]]

	cut_edges 		= [edge for edge in self.edges\
						if self.graph.nodes[edge[0]]["partition"] != self.graph.nodes[edge[1]]["partition"]]

	default_axes 	= plt.axes(frameon=True)
	

	if self.num_vertices > 20:
		stop 		= int(0.25 * self.num_vertices)
		pos			= nx.shell_layout(self.graph, [list(self.vertices)[:stop], list(self.vertices)[stop:]])
	else:
		pos 		= nx.circular_layout(self.graph)

	options 		= {"node_size": 1000,
						"width": 1.5,
						"pos": pos
						}



	nx.draw_networkx(self.graph, node_color=node_colors, edgelist=not_cut_edges, linewidths=1.5, edge_color= "black", font_size=14, alpha=1, ax=default_axes, **options)
	nx.draw_networkx_edges(self.graph, edgelist=cut_edges , edge_color= "brown",style='--', **options)

	labels 			= nx.get_edge_attributes(self.graph,'weight')
	# nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=labels)

	cut 			= 1.2
	xmax			= max(xx for xx,yy in pos.values())
	ymax			= max(yy for xx,yy in pos.values())
	xmin			= min(xx for xx,yy in pos.values())
	ymin			= min(yy for xx,yy in pos.values())

	xincrease 		= (cut - 1)*xmax
	yincrease 		= (cut - 1)*ymax

	plt.xlim(xmin - xincrease, cut*xmax)
	plt.ylim(ymin - yincrease, cut*ymax)
	
	figure_name 	= self.figure_path + "/"+ self.name +"_" + self.Params.Method + "_sol.png"

	plt.savefig(figure_name, dpi=300)
	plt.clf()


#================================================================================================
# Plot the decomposition tree
#================================================================================================
def plot_decomposition_tree(self):

	default_axes 	= plt.axes(frameon=True)
	position 		= nx.nx_pydot.pydot_layout(self.decomposition_tree, prog="dot")
	pos				= {key:(value[0], value[1]) for (key, value) in position.items()}

	options 		= {"font_size": 5,
						"node_color": "white",
						"edgecolors": "black",
						"linewidths": .6,
						"width": .6,
						}



	nx.draw_networkx(self.decomposition_tree.to_undirected(), ax=default_axes, pos=pos, **options)

	figure_name 	= self.figure_path + "/"+ self.name +"_" + self.Params.Method + "_decm_tree.png"

	cut 			= 1.1 
	xmax			= max(xx for xx,yy in pos.values())
	ymax			= max(yy for xx,yy in pos.values())
	xmin			= min(xx for xx,yy in pos.values())
	ymin			= min(yy for xx,yy in pos.values())

	xincrease 		= (cut - 1)*xmax
	yincrease 		= (cut - 1)*ymax

	plt.xlim(xmin - xincrease, cut*xmax)
	plt.ylim(ymin - yincrease, cut*ymax)

	plt.savefig(figure_name,  dpi=300)
	plt.clf()



#================================================================================================
# Plot the decomposition tree
#================================================================================================
def plot_qaoa_level_one(self, gammas, betas, objectives):

	objectives_min, objectives_max 	= -np.abs(objectives).max(), np.abs(objectives).max()

	fig, ax 						= plt.subplots()

	colormesh						= ax.pcolormesh(betas, gammas, objectives, cmap='RdBu', vmin=objectives_min, vmax=objectives_max, shading='auto')
	ax.set_title('QAOA$_1$-' + self.Params.Method)

	ax.axis([betas.min(), betas.max(), gammas.min(), gammas.max()])
	fig.colorbar(colormesh, ax=ax)

	figure_name 	= self.figure_path + "/"+ self.name +"_" + self.Params.Method + "_qaoa_level_one_heat_map.png"


	plt.ylabel(r'$\gamma/\pi$')
	plt.xlabel(r'$\beta/\pi$')


	plt.savefig(figure_name,  dpi=300)

	plt.clf()
	ax_3d 				= plt.axes(projection='3d')
	ax_3d.plot_surface(betas, gammas, objectives, rstride=1, cstride=1,
	                cmap='RdBu', edgecolor='none')
	ax_3d.set_title('QAOA$_1$-' + self.Params.Method)

	figure_name 	= self.figure_path + "/"+ self.name +"_" + self.Params.Method + "_qaoa_level_one_3d.png"

	plt.grid(False)

	# ax_3d.set_xticks([])
	# ax_3d.set_yticks([])
	# ax_3d.set_zticks([])


	ax_3d.xaxis.pane.fill = False
	ax_3d.yaxis.pane.fill = False
	ax_3d.zaxis.pane.fill = False


	ax_3d.xaxis.pane.set_edgecolor('w')
	ax_3d.yaxis.pane.set_edgecolor('w')
	ax_3d.zaxis.pane.set_edgecolor('w')

	ax_3d.set_xlabel(r'$\beta/\pi$', fontsize=10)
	ax_3d.set_ylabel(r'$\gamma/\pi$', fontsize=10)
	
	ax_3d.zaxis.set_rotate_label(False) 
	if self.Params.Method == "R-QUBO":
		ax_3d.set_zlabel(r'$\bar{q}$', fontsize=10, rotation = 0)
	elif self.Params.Method == "QUBO":
		ax_3d.set_zlabel(r'$q$', fontsize=10, rotation = 0)
	plt.savefig(figure_name,  dpi=300)
	plt.clf()

