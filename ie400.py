import numpy as np
import gurobipy as gp
from gurobipy import GRB
import time

# get the start time
st = time.time()

m = gp.Model("IE400")
PART = 'd'

stackPairItem = (("A","U"),("C","G"),("G","C"),("G","U"),("U","G"),("U","A"))
stackPairKey = [(a, b) for a in stackPairItem for b in stackPairItem]
    
# PARAMETERS
# S = "UGAUGGGUAUAAGACGAAGUUCGCCCAGUUGGCUCGAUUUGGUUGGUUGGCAGCUUACUACCUGGUUUCC"
#S = "UGAUGCGGUCAUC"
#S = "ACCCCUG"
S = "ACAAAAGU"
N = len(S)

matching_pairs = (('A', 'U'), ('C', 'G'), ('G', 'C'), ('U', 'A'))

energy_levels = {
    ('A', 'U'): -1.33,
    ('U', 'A'): -1.33, 
    ('G', 'C'): -1.45,
    ('C', 'G'): -1.45
}

# for part d
E = [
    [-1.1, -2.1, -2.2, -0.6],
    [-2.1, -2.4, -3.3, -1.4],
    [-2.2, -3.3, -3.4, -1.5],
    [-0.6, -1.4, -1.5, -0.3]
]
stack_pair_energy = {(matching_pairs[i], matching_pairs[j]): E[i][j] for i in range(4) for j in range(4)}


# DECISION VARIABLES
# x[i, j]
x = m.addMVar((N, N), vtype=GRB.BINARY, name="x")

# y[i, j]
if PART == 'd':
    y = m.addMVar((N, N), vtype=GRB.BINARY, name="y")


# MODEL
# objective function:
if PART == 'a':
    m.setObjective(x.sum(), GRB.MAXIMIZE)
elif PART in ('b', 'c'):
    m.setObjective(gp.quicksum(x[i, j] * energy_levels.get((S[i], S[j]), 0) for i in range(N)
                                                                            for j in range(N)), GRB.MINIMIZE)
elif PART == 'd':
    m.setObjective(gp.quicksum(y[i, j] * stack_pair_energy.get(((S[i], S[j]), (S[i+1], S[j-1])), 0) for i in range(N-1)
                                                                                                    for j in range(N)), GRB.MINIMIZE)

# at most one pairing per row (for each i)
m.addConstr(x.sum(axis=1) <= 1, name="row")

# at most one pairing per column (for each j)
m.addConstr(x.sum(axis=0) <= 1, name="col")

# if nucleotides do not match, x must be 0
m.addConstrs(x[i, j] == 0 for i in range(N)
                          for j in range(i + 1, N)
                          if (S[i], S[j]) not in matching_pairs)

# no overlapping
# if x[i, j] + x[i_, j_] = 2  then  (y[i, j]=1 or y[i-1, j+1]=1) and (y[i', j']=1 or y[i-1, j+1]=1) 
m.addConstrs(x[i, j] + x[i_, j_] <= 1 for i in range(N)
                                      for i_ in range(i + 1, N)
                                      for j in range(i_ + 1, N)
                                      for j_ in range(j + 1, N))

if PART != 'c': # a, b, d
    # x = 0 for any pair with less than 4 distance
    m.addConstrs(x[i, j] == 0 for i in range(N) 
                              for j in range(min(i + 4, N)))
else: # c
    m.addConstrs(x[i, j] == 0 for i in range(N) 
                              for j in range(min(i + 7, N)))

# if y[i, j] = 1  then  (x[i, j] = 1 and x[i+1, j-1] = 1)
if PART == 'd':
    m.addConstrs(2*y[i, j] <= x[i, j] + x[i+1, j-1] for i in range(N-1) 
                                                    for j in range(N))
    

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