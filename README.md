# Branch-Cut-and-Price Algorithm (BPCA) for the Bin Packing with Conflicts (BPC) Problem.

The solution was implemented using Julia JuMP, to be used with the Gurobi solver.

# Features:
  
- Subset Row cuts
- Ryan and Foster branching scheme
- Martello L2 lower bound calculation
- First Fit Decreasing (FFD) heuristic considering conflicts
- Simple labelling price (paralellized)
