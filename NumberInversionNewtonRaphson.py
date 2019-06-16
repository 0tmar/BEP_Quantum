from cQASM import *
from cQASMCompoundGates import ncx
from AdderQFT import cSUB
from QFT import QFT, iQFT
from MultiplierQFT import MUL, MULSUB


def invert_power(n, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamesd=None):
    """a: contains original number;
       b: will contain inverted power;
       c: contains ancilla qubit"""

    qna = buildnames(n=n, qubitnames=qubitnamesa, defaultname="a")
    qnb = buildnames(n=n, qubitnames=qubitnamesb, defaultname="b")
    qnc = buildnames(n=1, qubitnames=qubitnamec, defaultname="c")
    qnd = buildnames(n=n-2, qubitnames=qubitnamesd, defaultname="d")

    gates = []

    gates += [Qgate('cx', qna[0], qnb[n-1])]
    gates += [Qgate('cx', qna[0], qnc[0])]
    for i in range(n-1):
        gates += [Qgate('x', qnc[0])]
        gates += [Qgate('toffoli', qna[i+1], qnc[0], qnb[n-(i+2)])]
        gates += [Qgate('x', qnc[0])]
        gates += [Qgate('cx', qnb[n-(i+2)], qnc[0])]
    gates += ncx(qna=qna, qnb=qnc, qnc=qnd, invert=True, add_comment=False)
    gates += [Qgate('x', qnc[0])]
    gates += [Qgate('display')]

    return Qsubroutine(name='powerinvert', gates=gates)


def MULSUB_QFT_No_Parity(n, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamed=None):
    """a: eigenvalue, which will not be touched (size n);
       b: inverted power x0, which will be transformed to the first estimate x1 (size 2*n);
       c: control whether the subtraction is performed (size 1);
       d: ancilla qubit (size 1);
       (a, b, c, d) -> (a, b*(2-ab), c, d)"""

    qna = buildnames(n=n, qubitnames=qubitnamesa, defaultname="a")
    qnb = buildnames(n=2*n, qubitnames=qubitnamesb, defaultname="b")
    qnc = buildnames(n=1, qubitnames=qubitnamec, defaultname="c")
    qnd = buildnames(n=1, qubitnames=qubitnamed, defaultname="d")

    subroutines = []

    for i in range(n):
        cnotsubroutine = Qsubroutine(name='cnot', gates=[Qgate('cx', qnb[n-1-i], qnc[0])])
        qftsubroutine = QFT(n=n+1-i, qubitnames=qnb[(n - 1 - i):(2 * n - 2 * i)])
        subsubroutine = cSUB(na=n+1-i, nb=n-i, qubitnamesa=qnb[(n-1-i):(2*n-2*i)], qubitnamesb=qna[i:], qubitnamec=qnc,
                             qubitnamed=qnd)
        subsubroutine.gates.insert(0, Qgate('#', ' cSUB on {},...,{} by {},...,{} controlled by {}'.format(qnb[n-1-i], qnb[2*n-2*i-1], qna[i], qna[-1], qnc[0])))
        iqftsubroutine = iQFT(n=n+1-i, qubitnames=qnb[(n - 1 - i):(2 * n - 2 * i)])
        iqftsubroutine.gates.append(Qgate('display'))
        iqftsubroutine.gates.append(Qgate('prepz', qnc[0]))
        subroutines += [cnotsubroutine, qftsubroutine, subsubroutine, iqftsubroutine]
    # subroutines += [Qsubroutine(name='x', gates=[Qgate('x', qnc[0])])]

    return subroutines


class NUMINVcircuit(Qfunction):

    def __init__(self, inp="0", order=1):

        if not isinstance(inp, str):
            raise TypeError("input must be of type string")
        if order > 1:
            raise ValueError("order must be either 0 or 1 at this point")

        name = "Cao Number Inversion Estimator"
        n = len(inp)
        if abs(int(n)) is not n:
            raise ValueError("n must be integer")
        qna = buildnames(n=n, qubitnames="a")
        qnb = buildnames(n=2*n, qubitnames="b")
        qnc = buildnames(n=1, qubitnames="c")
        if order is 0:
            qubits = 3*n + 1
            qnd = []
            qn = qna + qnb + qnc
        elif order is 1:
            qubits = 3*n + 2
            qnd = buildnames(n=1, qubitnames="d")
            qn = qna + qnb + qnc + qnd
        else:
            raise ValueError("order must be either 0 or 1")

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(n):
            if inp[i] == "1":
                initgates += [Qgate("x", qna[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        subroutines = []
        subroutines += [initsubroutine]

        powinvsubroutine = invert_power(n=n, qubitnamesa=qna, qubitnamesb=qnb[:n], qubitnamec=qnc, qubitnamesd=qnb[n:(2*n-2)])
        subroutines += [powinvsubroutine]

        if order >= 1:
            mulsubsubroutines = MULSUB_QFT_No_Parity(n=n, qubitnamesa=qna, qubitnamesb=qnb, qubitnamec=qnc, qubitnamed=qnd)
            subroutines += mulsubsubroutines

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines += [resultsubroutine]

        super().__init__(name=name, qubits=qubits, subroutines=subroutines)
