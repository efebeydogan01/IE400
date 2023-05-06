import numpy as np
import gurobipy as gp
from gurobipy import GRB
import time

# get the start time
st = time.time()

m = gp.Model("PartA")
    
# PARAMETERS
S = "UGAUGGGUAUAAGACGAAGUUCGCCCAGUUGGCUCGAUUUGGUUGGUUGGCAGCUUACUACCUGGUUUCC"
N = len(S)

matching_pairs = (('A', 'U'), ('C', 'G'), ('G', 'C'), ('U', 'A'))

energy_levels = {
    ('A', 'U'): -1.33,
    ('U', 'A'): -1.33, 
    ('G', 'C'): -1.45,
    ('C', 'G'): -1.45
}

# DECISION VARIABLES
# x[i, j]: 1 if (i, j) paired
x = m.addMVar((N, N), vtype=GRB.BINARY, name="x")

# MODEL
# objective function:
m.setObjective(x.sum(), GRB.MAXIMIZE)


# at most one pairing per row (for each i)
m.addConstr(x.sum(axis=1) <= 1, name="row")
# at most one pairing per column (for each j)
m.addConstr(x.sum(axis=0) <= 1, name="col")

# one nucleotide cannot be in two pairings
m.addConstrs(x[i, j] + x[j, i_] <= 1 for i in range(N)
                                     for j in range(i+1, N)
                                     for i_ in range(j+1, N))

# matching constraint
# if nucleotides do not match, x must be 0
m.addConstrs(x[i, j] == 0 for i in range(N)
                          for j in range(i + 1, N)
                          if (S[i], S[j]) not in matching_pairs)

# overlapping constraint
m.addConstrs(x[i, j] + x[i_, j_] <= 1 for i in range(N)
                                      for j in range(i + 1, N)
                                      for i_ in range(i + 1, j)
                                      for j_ in range(j + 1, N))
# distance constraint
# x = 0 for any pair with less than 4 distance
m.addConstrs(x[i, j] == 0 for i in range(N) 
                          for j in range(min(i + 4, N)))
    
# solve the problem
m.optimize()

# print solution
if m.status == GRB.OPTIMAL:
    print("Optimal solution found")
    pairs = []
    matches = []
    for i in range(N):
        for j in range(N):
            if x[i,j].X == 1:
                pairs.append((i,j))
                matches.append((S[i], S[j]))
    print("Pairs:", pairs)
    print("Matches:", matches)
    print("pairs match:", all(match in matching_pairs for match in matches))
    print("Number of pairs:", len(pairs))
    print("Optimal value:", m.ObjVal)
    print("Runtime:", m.Runtime)
else:
    print("No solution found")

# print the execution time
elapsed_time = time.time() - st
print('Execution time:', elapsed_time, 'seconds')