import numpy as np
import gurobipy as gp
from gurobipy import GRB
import time

# get the start time
st = time.time()

m = gp.Model("IE400f")
    
# PARAMETERS
# S = "UGAUGGGUAUAAGACGAAGUUCGCCCAGUUGGCUCGAUUUGGUUGGUUGGCAGCUUACUACCUGGUUUCC"
S = "UGAUGCGGUCAUC"
#S = "AACCUUUUCGGAA"
#S = "ACCCCUGGG"
#S = "ACAAAAGU"
#S = "AACCCUU"
N = len(S)

matching_pairs = (('A', 'U'), ('C', 'G'), ('G', 'C'), ('U', 'A'))

E = [
    [-1.1, -2.1, -2.2, -0.6],
    [-2.1, -2.4, -3.3, -1.4],
    [-2.2, -3.3, -3.4, -1.5],
    [-0.6, -1.4, -1.5, -0.3]
]
stack_pair_energy = {(matching_pairs[i], matching_pairs[j]): E[i][j] for i in range(4) for j in range(4)}


# DECISION VARIABLES
# e[i, j]
e = m.addMVar((N, N), lb=-GRB.INFINITY, ub=0, name="e")

a = m.addMVar((N, N, N), lb=-GRB.INFINITY, ub=0, name="a")
b = m.addMVar((N, N), lb=-GRB.INFINITY, ub=0, name="b")

# x[i, j]: 1 if (i, j) paired
x = m.addMVar((N, N), vtype=GRB.BINARY, name="x")

# MODEL
"""
e[i, j]: min total energy constructed from the substring from i to j

            | 0                                                                                 if i >= j - 4
e[i, j] =   | min_(i<=k<j){e[i, k] + e[k+1, j]}                                                 if (i, j) or (i+1, j-1) not matched
            | min(min_(i<=k<j){e[i, k] + e[k+1, j]}, SPE[(i, j), (i+1, j-1)] + e[i+2, j-2])     otherwise
"""
# objective function:
m.setObjective(e[0, N-1], GRB.MINIMIZE)

m.addConstrs(e[i, j] == 0 for i in range(N)
                          for j in range(min(i + 4, N)))

for i in range(N):
    for j in range(i+4, N):
        t1, t2 = (S[i], S[j]), (S[i+1], S[j-1])
        m.addConstrs(a[i, k, j] == e[i, k] + e[k+1, j] for k in range(i, j))

        if t1 in matching_pairs and t2 in matching_pairs:
            m.addConstr(b[i, j] == stack_pair_energy[(t1, t2)] + e[i+1, j-1])
            m.addConstr(e[i, j] == gp.min_([a[i, k, j] for k in range(i, j)] + [b[i, j]], 0))
        else:
            m.addConstr(e[i, j] == gp.min_([a[i, k, j] for k in range(i, j)], 0))

# solve the problem
m.optimize()

# print solution
if m.status == GRB.OPTIMAL:
    print("Optimal solution found")
    pairs = []
    matches = []
    print(e.X)
    for i in range(N):
        for j in range(N):
            if e[i,j].X == 1:
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
