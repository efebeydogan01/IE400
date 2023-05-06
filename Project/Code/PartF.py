import numpy as np
import gurobipy as gp
from gurobipy import GRB
import time

# get the start time
st = time.time()
    
# PARAMETERS
S = "UGAUGGGUAUAAGACGAAGUUCGCCCAGUUGGCUCGAUUUGGUUGGUUGGCAGCUUACUACCUGGUUUCC"
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
e = np.zeros((N, N))
b = np.zeros((N, N))

# MODEL
"""
e[i, j]: min total energy constructed from the substring from i to j

            | 0                                                     if i > j - 4
e[i, j] =   | min(min_(i<=k<j){e[i, k] + e[k+1, j]}, b[i, j])       if (i, j) matching
            | min_(i<=k<j){e[i, k] + e[k+1, j]}                     otherwise                        


            
b[i, j]: min total energy constructed from the substring from i to j if (i, j) is paired

            | 0                                                     if i > j - 4  or  (i, j) not matching
b[i, j] =   | SPE[(i, j), (i+1, j-1)] + b[i+1, j-1]                 if j >= i + 6  and (i+1, j-1) matching
            | e[i+1, j-1]                                           otherwise
"""

for j in range(4, N):
    for i in range(j-4, -1, -1):
        pair, next_pair = (S[i], S[j]), (S[i+1], S[j-1])

        if pair not in matching_pairs:
            b[i, j] = 0
        elif j >= i + 6 and next_pair in matching_pairs:
            b[i, j] = stack_pair_energy[(pair, next_pair)] + b[i+1, j-1]
        else:
            b[i, j] = e[i+1, j-1]

        if pair in matching_pairs:
            e[i, j] = min(b[i, j], min(e[i, k] + e[k+1, j] for k in range(i, j)))
        else:
            e[i, j] = min(e[i, k] + e[k+1, j] for k in range(i, j))

# get results
pairs = []
stacks = []

# for backtracking
def get_pairs(i, j):
    global pairs, stacks
    if e[i,j] == 0:
        return
    if e[i,j] == b[i,j]:
        pairs.append((i,j))

        while b[i, j] != e[i+1, j-1] and b[i, j] != 0:
            pair_tuple = ((S[i], S[j]), (S[i+1], S[j-1]))
            pairs.append((i+1,j-1))
            stacks.append((pair_tuple, stack_pair_energy[pair_tuple]))
            i = i + 1
            j = j - 1
        if b[i, j] == e[i+1, j-1]:
            get_pairs(i+1, j-1)
    else:
        k = np.argmin([e[i, k] + e[k+1, j] for k in range(i, j)])
        get_pairs(i, i+k)
        get_pairs(i+k+1, j)

get_pairs(0, N-1)

matches = [(S[i], S[j]) for i, j in pairs]

print("Pairs:", pairs)
print("Matches:", matches)
print("Stacks:", stacks)
print("Total Energy:", sum(s[1] for s in stacks))
print("pairs match:", all(match in matching_pairs for match in matches))
print("Number of pairs:", len(pairs))
print("Optimal value:", round(e[0, N-1], 2))

# print the execution time
elapsed_time = time.time() - st
print('Execution time:', elapsed_time, 'seconds')