import itertools
import math
import numpy as np
from scipy.special import poch  # The Pochhammer symbol (x)_n = Gamma(x+n)/Gamma(x) = x*(x+1)*...*(x+n-1)
from cQASM import *
from cQASMCompoundGates import cRy, ccRy


def x_to_the_p_rot_combinations_and_angles(n, p):
    """Returns the controlling qubit combinations and corresponding Ry angles for the n-unique-controlling qubits in
    the Ry(c*x^p) operation"""

    a = list(range(n))
    ctrl_rot_array = [[]]*min(n, p)
    for l in range(1, min(n, p) + 1):
        # Combinations of qubits to control the l-fold-controlled-Ry gates
        combinations = list(set(itertools.combinations(a, l)))

        # Possible internal permutations of "how many times" each qubit controls the Ry gate
        # (p = boxes, l = balls)
        rng = list(range(1, p + 1)) * l
        permutations = list(set(i for i in itertools.permutations(rng, l) if sum(i) == p))
        permutations_with_relative_frequencies = list([perm, int(math.factorial(p) / (np.prod(list(math.factorial(val) for val in perm))))] for perm in permutations)

        ctrl_rot_array_temp = [([0], 0)]*len(combinations)
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


def p_fold_controlled_Ry(qnc, qnr, qna, angle):
    """Returns the gates for performing a p-fold-controlled-Ry(theta) operation. Here, qnc are the controlling qubits,
    qnr is the qubit being rotated, and qna are the possible ancillae, of which at least max(0, p-2) are required"""

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


def c_x_to_the_p_rot_gates(qnc, qnr, qna, c, p):
    """Returns the gates for performing an Ry(c*x^p) operation, where x is stored in qnc, qnr is rotated, and qna are
    ancillae, of which at least max(0,p-2) are required"""

    n = len(qnc)
    ctrl_rot_array = x_to_the_p_rot_combinations_and_angles(n, p)

    gates = []
    gates += [Qgate("#", "#### Performs Ry(c*(x^p)) for c={} and p={}".format(c, p))]
    gates += [Qgate()]
    for i in range(min(len(qnc), p)):
        gates += [Qgate("#", "## Rotations controlled by {} qubit(s)".format(i+1))]
        gates += [Qgate()]
        for j in range(len(ctrl_rot_array[i])):
            comb = ctrl_rot_array[i][j][0]
            theta = ctrl_rot_array[i][j][1]
            qnc_temp = [qnc[i] for i in comb]
            gates += [Qgate("#", " Rotation c*theta with c = {} and theta = {}, controlled by {}".format(c, theta, comb))]
            gates += p_fold_controlled_Ry(qnc=qnc_temp, qnr=qnr, qna=qna, angle=2 * c * theta)
            gates += [Qgate()]
    return gates


def arcsin_taylor_factor(p):
    """Returns the (1+2p)th Taylor expansion factor of arcsin(x) around x=0, which is a_p = (1/2)_p / (1+2p)*p!"""

    return poch(1 / 2, p) / ((1 + 2 * p) * math.factorial(p))


def ancilla_rotation_subroutines(qnc, qnr, qna, c=1, k=2):
    """Returns the subroutines for performing a rotation from |0> to [(c*x)|1> + sqrt(1-(c*x)^2)|0>], by approximating
    an arcsin up to the m-th term"""

    n = len(qnc)
    subroutines = []

    for i in range(k + 1):
        a_i = arcsin_taylor_factor(i)
        c_i = c**(1+2*i)
        subroutines += [Ry_c_x_to_the_p(n=n, c=a_i * c_i, p=1 + 2 * i, qubitnamesc=qnc, qubitnamer=qnr, qubitnamesa=qna)]

    return subroutines


class Ry_c_x_to_the_p(Qsubroutine):
    """Quantum subroutine for performing an Ry(c*x^p) operation, where x is an n qubit number"""

    def __init__(self, n=4, c=1, p=3, qubitnamesc=None, qubitnamer=None, qubitnamesa=None):

        qnc = buildnames(n, qubitnamesc, "c")
        qnr = buildnames(1, qubitnamer, "r")
        if qubitnamesa is not None:
            if len(qubitnamesa) >= max(0, min(n, p) - 2):
                qna = buildnames(len(qubitnamesa), qubitnamesa, "a")
            else:
                raise IndexError("Input for qubitnamesa is too short: {} instead of at least {}".format(len(qubitnamesa), max(0, min(n, p) - 2)))
        else:
            qna = buildnames(max(0, min(n, p) - 2), qubitnamesa, "a")
        self.qubitnamesc = qnc
        self.qubitnamer = qnr
        self.qubitnamesa = qna

        gates = c_x_to_the_p_rot_gates(qnc=qnc, qnr=qnr, qna=qna, c=c, p=p)

        super().__init__(name="ry_c_x_to_the_p", gates=gates)


class Ry_c_x_to_the_p_circuit(Qfunction):
    """Quantum circuit for performing an Ry(c*x^p) operation"""

    def __init__(self, inp="0", c=1, p=3):
        name = "Ry(c*(x^p)) rotation for c={} and p={}".format(c, p)
        n = len(inp)
        qubits = n + max(0, min(n, p) - 2) + 1
        rysubroutine = Ry_c_x_to_the_p(n=n, c=c, p=p)
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


class AncillaRotation(Qsubroutine):
    """Quantum subroutine for performing the Ancilla Rotation subroutine in the HHL Quantum Linear Solver Algorithm"""

    def __init__(self, n=4, c=1, k=3, qubitnamesc=None, qubitnamer=None, qubitnamesa=None):

        qnc = buildnames(n, qubitnamesc, "c")
        qnr = buildnames(1, qubitnamer, "r")
        qna = buildnames(max(0, min(n, 1 + 2 * k) - 2), qubitnamesa, "a")
        self.qubitnamesc = qnc
        self.qubitnamer = qnr
        self.qubitnamesa = qna

        gates = []
        gates += [Qgate("#", "###### Ancilla rotation using the first m terms of arcsin(x) and c={}".format(c))]
        gates += [Qgate()]
        gates += [i for j in ancilla_rotation_subroutines(qnc=qnc, qnr=qnr, qna=qna, c=c, k=k) for i in j.gates]

        super().__init__(name="ancilla_rotation", gates=gates)


class AncillaRotationCircuit(Qfunction):
    """Quantum circuit for performing the Ancilla Rotation subroutine in the HHL Quantum Linear Solver Algorithm"""

    def __init__(self, inp="0", c=1, k=3):
        name = "Ancilla rotation for c={}, approximated using the m={} first terms of the taylor expansion of arcsin(x)".format(c, k)
        n = len(inp)
        qubits = n + max(0, min(n, 1 + 2 * k) - 2) + 1
        ancilla_rotation_subroutine = AncillaRotation(n=n, c=c, k=k)
        qnc = ancilla_rotation_subroutine.qubitnamesc
        qnr = ancilla_rotation_subroutine.qubitnamer
        qna = ancilla_rotation_subroutine.qubitnamesa
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

        ancrotsubroutines = ancilla_rotation_subroutines(qnc=qnc, qnr=qnr, qna=qna, c=c, k=k)

        resultgates = [Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = []
        subroutines += [initsubroutine]
        subroutines += ancrotsubroutines
        subroutines += [resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)
