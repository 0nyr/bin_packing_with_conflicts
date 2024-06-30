# Branch-Cut-and-Price Algorithm (BPCA) for the Bin Packing with Conflicts (BPC) Problem.

The solution was implemented using Julia JuMP and Gurobi.

# Features:
  
- Subset Row cuts
- Ryan and Foster branching scheme
- Martello L2 lower bound calculation
- First Fit Decreasing (FFD) heuristic considering conflicts
- Pricing heuristic (relaxed LP pricing, rounded to integer)
- Simple labelling price (paralellized)
