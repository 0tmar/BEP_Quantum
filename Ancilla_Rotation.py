import itertools
import math
import numpy as np


n = 4
k = 3

a = list(range(n))
for l in range(1, k+1):
    # Combinations of qubits to control the l-fold-controlled-Ry gates
    combinations = set(itertools.combinations(a, l))

    # Possible internal permutations of "how many times" each qubit controls the Ry gate
    # (k = boxes, l = balls)
    rng = list(range(1, k+1)) * l
    permutations = list(set(i for i in itertools.permutations(rng, l) if sum(i) == k))
    permutations_with_relative_frequencies = list([perm, int(math.factorial(k) / (np.prod(list(math.factorial(val) for val in perm))))] for perm in permutations)

    print(permutations_with_relative_frequencies)

    for comb in combinations:
        comb = np.array(comb)
        angle = 0
        for i in range(len(permutations)):
            perm_plus_freq = permutations_with_relative_frequencies[i]
            perm = np.array(perm_plus_freq[0])
            rel_freq = perm_plus_freq[1]

            angle += rel_freq*np.sum(2.**(-comb*perm))
        print(comb)
        print(angle)

    print("")


