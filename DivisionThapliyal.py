from cQASM import *
from AdderMunozCoreas import *


class DIV(Qsubroutine):

    def __init__(self, n=1, m=None, qubitnamesn=None, qubitnameso=None, qubitnamesd=None, sign=1):

        if not (sign == -1 or sign == 1):
            raise ValueError("Variable 'sign' should be either 1 or -1, not {}".format(sign))

        def iteration(n, qny, qnz, qnd, i, sign=1):

            gates = []

            if sign == 1:

                subsubroutine = SUB(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamez=qnz, do_overflow=True, do_ctrl=False, subtype="b-a")
                caddsubroutine = ADD(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamectrl=qnz, do_overflow=False, do_ctrl=True)

                gates += [Qgate("#", "## Iteration {}, where Y = [{}, ..., {}] and Z = {}".format(i, qny[0], qny[-1], qnz))]
                gates += [Qgate()]
                gates += [Qgate("#", " SUB: Y - D, on Y")]
                gates += subsubroutine.gates
                # gates += [Qgate()]
                # gates += [Qgate("#", " cx on Z, ctrl'd by Y[-1]")]
                # gates += [Qgate("cx", qny[-1], qnz)]
                gates += [Qgate()]
                gates += [Qgate("#", " cADD: Y (+ D), on Y, ctrl'd by Z")]
                gates += caddsubroutine.gates
                gates += [Qgate()]
                gates += [Qgate("#", " x on Z")]
                gates += [Qgate("x", qnz)]

            else:

                unsubsubroutine = ADD(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamez=qnz, do_overflow=True, do_ctrl=False)
                uncaddsubroutine = SUB(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamectrl=qnz, do_overflow=False, do_ctrl=True, subtype="b-a")

                gates += [Qgate("#", "## Un-iteration {}, where Y = [{}, ..., {}] and R = {}".format(i, qny[0], qny[-1], qnz))]
                gates += [Qgate()]
                gates += [Qgate("#", " x on Z")]
                gates += [Qgate("x", qnz)]
                gates += [Qgate()]
                gates += [Qgate("#", " un-cADD: Y (+ D), on Y, ctrl'd by Z")]
                gates += uncaddsubroutine.gates
                gates += [Qgate()]
                # gates += [Qgate()]
                # gates += [Qgate("#", " cx on Z, ctrl'd by Y[-1]")]
                # gates += [Qgate("cx", qny[-1], qnz)]
                gates += [Qgate("#", " un-SUB: Y - D, on Y")]
                gates += unsubsubroutine.gates

            return gates

        if m is None:
            m = n
        elif m > n:
            raise ValueError("Value of m must be smaller than or equal to n")

        qnn = buildnames(n, qubitnamesn, "n")
        qno = buildnames(m, qubitnameso, "o")
        qnd = buildnames(m, qubitnamesd, "d")

        self.qubitnamesn = qnn
        self.qubitnameso = qno
        self.qubitnamesd = qnd

        gates = []

        if sign == 1:
            gates += [Qgate("#", "#### Thapliyal division with N = [{}, ..., {}], O = [{}, ..., {}] and D = [{}, ..., {}]". format(qnn[0], qnn[-1], qno[0], qno[-1], qnd[0], qnd[-1]))]
            rng = list(range(n))
        else:
            gates += [Qgate("#", "#### un-Thapliyal division with N = [{}, ..., {}], O = [{}, ..., {}] and D = [{}, ..., {}]". format(qnn[0], qnn[-1], qno[0], qno[-1], qnd[0], qnd[-1]))]
            rng = list(reversed(range(n)))

        gates += [Qgate()]

        for i in rng:

            if i < m:
                qny = qnn[-i-1:] + qno[:m-i-1]
                qnztemp = qno[m-i-1]
            else:
                qny = qnn[-i-1:m-i-1]
                qnztemp = qnn[m-i-1]

            gates += iteration(n=m, qny=qny, qnd=qnd, qnz=qnztemp, i=i+1, sign=sign)

            if i is not rng[-1]:
                gates += [Qgate()]

        super().__init__(name="thapliyal_division", gates=gates)


class DIVcircuit(Qfunction):

    def __init__(self, inp_n="0", inp_d="0"):
        name = "Thapliyal Division"
        nn = len(inp_n)
        nd = len(inp_d)
        n = max(nn, nd)
        inp_n = (n-nn)*"0" + inp_n
        inp_d = (n-nd)*"0" + inp_d
        qubits = 3*n
        divsubroutine = DIV(n=n)
        qnn = divsubroutine.qubitnamesn
        qno = divsubroutine.qubitnameso
        qnd = divsubroutine.qubitnamesd
        qn = qnn + qno + qnd

        if not isinstance(inp_n, str):
            raise TypeError("input for inp_n must be of type string")
        if not isinstance(inp_d, str):
            raise TypeError("input for inp_d must be of type string")

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(n):
            if inp_n[-i - 1] == "1":
                initgates += [Qgate("x", qnn[i])]
        for i in range(n):
            if inp_d[-i - 1] == "1":
                initgates += [Qgate("x", qnd[i])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, divsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


class DivUnequalCircuit(Qfunction):

    def __init__(self, inp_n="0", inp_d="0"):
        if not isinstance(inp_n, str):
            raise TypeError("input for inp_n must be of type string")
        if not isinstance(inp_d, str):
            raise TypeError("input for inp_d must be of type string")

        name = "Thapliyal Division"
        nn = len(inp_n)
        no = len(inp_d)
        n = nn
        m = no
        if inp_d[0] == '1' and n > m:
            inp_d = "0" + inp_d
            m += 1
        qubits = n + 2*m
        divsubroutine = DIV(n=n, m=m)
        qnn = divsubroutine.qubitnamesn
        qno = divsubroutine.qubitnameso
        qnd = divsubroutine.qubitnamesd
        qn = qnn + qno + qnd

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(n):
            if inp_n[-i-1] == "1":
                initgates += [Qgate("x", qnn[i])]
        for i in range(m):
            if inp_d[-i-1] == "1":
                initgates += [Qgate("x", qnd[i])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, divsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


if __name__ == "__main__":

    pass
