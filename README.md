# The max k-cut problem on classical and quantum solvers

Software package for the papers "[The max k-cut problem on classical and quantum solvers](https://www.sciencedirect.com/science/article/abs/pii/S0167637723001293)" and "A folding preprocess for the max k-cut problem" by Ramin Fakhimi, Hamidreza Validi, Illya V. Hicks, Tamas Terlaky, and Luis F. Zuluaga.

The max k-cut problem is easy to be explained, but hard to be solved: How to color the nodes of a graph with at most k colors such that the number of edges with distinct colored endpoints is maximized?

![Figure 1](readme_images/input_graph.png?raw=true "Input graph")
![Figure 2](readme_images/solution_max_3-cut.png?raw=true "An optimal solution for the max 3-cut problem")

We study four classical mixed integer linear optimization models of the max k-cut problem and provide theoretical and computational comparisons between them. As the classical models cannot be fed to quantum machines directly, we propose two quadratic unconstrained binary optimization (QUBO) models with tight penalty coefficients. 


## Install

1- Use the following comands to install the max_k_cut package

```
git clone https://github.com/QCOL-LU/MaxKcut.git
python3.9 -m venv env
source env/bin/activate
pip3.9 install -e MaxKcut/max_k_cut
pip3 install qiskit numpy scipy matplotlib networkx qiskit-ibm-runtime
```

2- Follow the instruction provided [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-) to install Gurobi solver.

3- Follow the instruction provided [here](https://docs.mosek.com/latest/install/installation.html) to install Mosek solver.

## Run

1- Import ```max_k_cut``` to the code

2- Generate a networkx graph for the problem

3- Generate an instance object by ```my_instance = Instance(graph, name_specifier=name)```

4- Set the parameters of the algorithms as follows (all the paramteres are avilable in ```main\ParametersDefault.py```)

```my_instance.Params.Num_Partitions = num_partitions```

5- Execute ```my_instance.solve()```

