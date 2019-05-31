import itertools
import math
import numpy as np
from cQASM import *
from cQASM_Compound_Gates import cRy, ccRy


def x_to_the_k_rot_combinations_and_angles(n, k):
    a = list(range(n))
    ctrl_rot_array = [None]*min(n, k)
    for l in range(1, min(n, k)+1):
        # Combinations of qubits to control the l-fold-controlled-Ry gates
        combinations = list(set(itertools.combinations(a, l)))

        # Possible internal permutations of "how many times" each qubit controls the Ry gate
        # (k = boxes, l = balls)
        rng = list(range(1, k+1)) * l
        permutations = list(set(i for i in itertools.permutations(rng, l) if sum(i) == k))
        permutations_with_relative_frequencies = list([perm, int(math.factorial(k) / (np.prod(list(math.factorial(val) for val in perm))))] for perm in permutations)

        ctrl_rot_array_temp = [None]*len(combinations)
        for j in range(len(combinations)):
            comb = np.array(combinations[j])
            angle = 0
            for i in range(len(permutations)):
                perm_plus_freq = permutations_with_relative_frequencies[i]
                perm = np.array(perm_plus_freq[0])
                rel_freq = perm_plus_freq[1]
                angle += rel_freq*np.prod(2.**(-comb*perm))
            ctrl_rot_array_temp[j] = (list(combinations[j]), angle)
        ctrl_rot_array[l-1] = ctrl_rot_array_temp
    return ctrl_rot_array


def k_fold_controlled_Ry(qnc, qnr, qna, angle):
    if (len(qnc) >= 3) and (len(qnc) > len(qna)+2):
        raise IndexError("The size of qna must be at least as large as len(qnc)-2 for len(qnc)>=3")
    nc = len(qnc)
    gates = []
    if nc == 0:
        raise IndexError("Size of qnc must at least be one")
    elif nc == 1:
        gates += cRy(qna=qnc[0], qnb=qnr[0], theta=angle, add_comment=False)
    elif nc == 2:
        gates += ccRy(qna=qnc[0], qnb=qnc[1], qnc=qnr[0], theta=angle, add_comment=False)
    else:
        gates += [Qgate("toffoli", qnc[0], qnc[1], qna[0])]
        for i in range(nc-3):
            gates += [Qgate("toffoli", qnc[i+2], qna[i], qna[i+1])]
        gates += ccRy(qna=qnc[-1], qnb=qna[nc-3], qnc=qnr[0], theta=angle, add_comment=False)
        for i in reversed(range(nc-3)):
            gates += [Qgate("toffoli", qnc[i+2], qna[i], qna[i+1])]
        gates += [Qgate("toffoli", qnc[0], qnc[1], qna[0])]
    return gates


def c_x_to_the_k_rot_gates(qnc, qnr, qna, c, k):
    n = len(qnc)
    ctrl_rot_array = x_to_the_k_rot_combinations_and_angles(n, k)
    c3 = c**3

    gates = []
    for i in range(k):
        gates += [Qgate("#", "## Rotations controlled by {} qubit(s)".format(i))]
        gates += [Qgate()]
        for j in range(len(ctrl_rot_array[i])):
            comb = ctrl_rot_array[i][j][0]
            theta = ctrl_rot_array[i][j][1]
            qnc_temp = [qnc[i] for i in comb]
            gates += [Qgate("#", " Rotation c*theta with c = {} and theta = {}, controlled by {}".format(c, theta, comb))]
            gates += k_fold_controlled_Ry(qnc=qnc_temp, qnr=qnr, qna=qna, angle=2*c3*theta)
            gates += [Qgate()]
    return gates


class Ry_cx_to_the_k(Qsubroutine):

    def __init__(self, n=4, c=1, k=3, qubitnamesc=None, qubitnamer=None, qubitnamesa=None):

        qnc = buildnames(n, qubitnamesc, "c")
        qnr = buildnames(1, qubitnamer, "r")
        qna = buildnames(max(0, n-2), qubitnamesa, "a")
        self.qubitnamesc = qnc
        self.qubitnamer = qnr
        self.qubitnamesa = qna

        gates = c_x_to_the_k_rot_gates(qnc=qnc, qnr=qnr, qna=qna, c=c, k=k)

        super().__init__(name="ry_cx_to_the_k", gates=gates)


class Ry_cx_to_the_k_circuit(Qfunction):

    def __init__(self, inp="0", c=1, k=3):
        name = "Ry(c*(x^k) rotation for c={}, k={})".format(c, k)
        n = len(inp)
        qubits = n + max(0, n-2) + 1
        rysubroutine = Ry_cx_to_the_k(n=n, c=c, k=k)
        qnc = rysubroutine.qubitnamesc
        qnr = rysubroutine.qubitnamer
        qna = rysubroutine.qubitnamesa
        qn = []
        qn += qnc
        qn += qna
        qn += qnr

        if not isinstance(inp, str):
            raise TypeError("input must be of type string")

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(n):
            if inp[i] == "1":
                initgates += [Qgate("x", qnc[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        resultgates = [Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, rysubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)
