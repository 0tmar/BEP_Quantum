from cQASM import *
from AdderMunozCoreas import *


class DIV(Qfunction):

    def __init__(self, n=1, m=None, qubitnamesq=None, qubitnamesr=None, qubitnamesd=None, qubitnamec=None):

        def iteration(n, qny, qnr, qnd, i):

            subsubroutine = SUB(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamez=qnr, do_overflow=True, do_ctrl=False, subtype="b-a")
            caddsubroutine = ADD(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamectrl=qnr, do_overflow=False, do_ctrl=True)

            gates = []
            gates += [Qgate("#", "## Iteration {}, where Y = [{}, ..., {}] and R = {}".format(i, qny[0], qny[-1], qnr))]
            gates += [Qgate()]
            gates += [Qgate("#", " SUB: Y - D, on Y")]
            gates += subsubroutine.gates
            # gates += [Qgate()]
            # gates += [Qgate("#", " cx on R, ctrl'd by Y[-1]")]
            # gates += [Qgate("cx", qny[-1], qnr)]
            gates += [Qgate()]
            gates += [Qgate("#", " cADD: Y (+ D), on Y, ctrl'd by R")]
            gates += caddsubroutine.gates
            gates += [Qgate()]
            gates += [Qgate("#", " x on R")]
            gates += [Qgate("x", qnr)]

            return Qsubroutine(name="iteration", gates=gates)

        if m is None:
            m = n
        elif m > n:
            raise ValueError("Value of m must be smaller than or equal to n")

        qnq = buildnames(n, qubitnamesq, "q")
        qnr = buildnames(m, qubitnamesr, "r")
        qnd = buildnames(m, qubitnamesd, "d")

        self.qubitnamesq = qnq
        self.qubitnamesr = qnr
        self.qubitnamesd = qnd

        subroutines = []

        for i in range(n):

            if i < m:
                qny = qnq[-i-1:] + qnr[:m-i-1]
                qnrtemp = qnr[m-i-1]
            else:
                qny = qnq[-i-1:m-i-1]
                qnrtemp = qnq[m-i-1]
            subroutines += [iteration(n=m, qny=qny, qnd=qnd, qnr=qnrtemp, i=i+1)]

        super().__init__(name="Thapliyal Division", subroutines=subroutines)


class DIVcircuit(Qfunction):

    def __init__(self, inp_n="0", inp_d="0"):
        name = "Thapliyal Division"
        na = len(inp_n)
        nb = len(inp_d)
        n = max(na, nb)
        inp_n = (n-na)*"0" + inp_n
        inp_d = (n-nb)*"0" + inp_d
        qubits = 3*n
        divfunction = DIV(n=n)
        qnq = divfunction.qubitnamesq
        qnr = divfunction.qubitnamesr
        qnd = divfunction.qubitnamesd
        qn = qnq + qnr + qnd

        if not isinstance(inp_n, str):
            raise TypeError("input for inp_n must be of type string")
        if not isinstance(inp_d, str):
            raise TypeError("input for inp_d must be of type string")

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(n):
            if inp_n[-i - 1] == "1":
                initgates += [Qgate("x", qnq[i])]
        for i in range(n):
            if inp_d[-i - 1] == "1":
                initgates += [Qgate("x", qnd[i])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        divsubroutines = divfunction.subroutines

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine] + divsubroutines + [resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


class DivUnequalCircuit(Qfunction):

    def __init__(self, inp_n="0", inp_d="0"):
        if not isinstance(inp_n, str):
            raise TypeError("input for inp_n must be of type string")
        if not isinstance(inp_d, str):
            raise TypeError("input for inp_d must be of type string")

        name = "Thapliyal Division"
        na = len(inp_n)
        nb = len(inp_d)
        n = na
        m = nb
        if inp_d[0] == '1' and n > m:
            inp_d = "0" + inp_d
            m += 1
        qubits = n + 2*m
        divfunction = DIV(n=n, m=m)
        qnq = divfunction.qubitnamesq
        qnr = divfunction.qubitnamesr
        qnd = divfunction.qubitnamesd
        qn = qnq + qnr + qnd

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(n):
            if inp_n[-i-1] == "1":
                initgates += [Qgate("x", qnq[i])]
        for i in range(m):
            if inp_d[-i-1] == "1":
                initgates += [Qgate("x", qnd[i])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        divsubroutines = divfunction.subroutines

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine] + divsubroutines + [resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


if __name__ == "__main__":

    pass
