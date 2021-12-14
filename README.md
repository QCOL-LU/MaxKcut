# The max k-cut problem on classical and quantum solvers

Software package for the [paper](https://engineering.lehigh.edu/sites/engineering.lehigh.edu/files/_DEPARTMENTS/ise/pdf/tech-papers/21/21T_007.pdf) "The max k-cut problem on classical and quantum solvers" by Ramin Fakhimi, Hamidreza Validi, Illya V. Hicks, Tamas Terlaky, and Luis F. Zuluaga.

The max k-cut problem is easy to be explained, but hard to be solved: How to color the nodes of a graph with at most k colors such that the number of edges with distinct colored endpoints is maximized?

![Figure 1](input_graph.png?raw=true "Input graph")
![Figure 2](solution_max_3-cut.png?raw=true "An optimal solution for the max 3-cut problem")

We study four classical mixed integer linear optimization models of the max k-cut problem and provide theoretical and computational comparisons between them. As the classical models cannot be fed to quantum machines directly, we propose two quadratic unconstrained binary optimization (QUBO) models with tight penalty coefficients. 
