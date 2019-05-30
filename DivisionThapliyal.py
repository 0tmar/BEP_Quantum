from cQASM import *
from AdderCuccaro import *


class DIV(Qfunction):

    def __init__(self, n=1, qubitnamesq=None, qubitnamesr=None, qubitnamesd=None, qubitnamec=None):

        def iteration(n, qny, qnr, qnd, qnc, i):

            subsubroutine = SUB(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamec=qnc, qubitnamez=qnr, do_overflow=True, subtype="b-a")
            caddsubroutine = cADD(n=n, qubitnamesa=qnd, qubitnamesb=qny, qubitnamec=qnc, qubitnamectrl=qnr, do_overflow=False)

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

        qnq = buildnames(n, qubitnamesq, "q")
        qnr = buildnames(n, qubitnamesr, "r")
        qnd = buildnames(n, qubitnamesd, "d")
        qnc = buildnames(1, qubitnamec, "c")
        qnc = qnc[0]

        self.qubitnamesq = qnq
        self.qubitnamesr = qnr
        self.qubitnamesd = qnd
        self.qubitnamec = qnc

        subroutines = []

        for i in range(n):

            qny = qnq[-i-1:] + qnr[:n-i-1]
            qnrtemp = qnr[n-i-1]
            subroutines += [iteration(n=n, qny=qny, qnd=qnd, qnr=qnrtemp, qnc=qnc, i=i+1)]

        super().__init__(name="Thapliyal Division", subroutines=subroutines)


class DIVcircuit(Qfunction):

    def __init__(self, inp_n="0", inp_d="0"):
        name = "Thapliyal Division"
        na = len(inp_n)
        nb = len(inp_d)
        n = max(na, nb)
        inp_n = (n - na) * "0" + inp_n
        inp_d = (n - nb) * "0" + inp_d
        qubits = 3*n + 1
        divfunction = DIV(n=n)
        qnq = divfunction.qubitnamesq
        qnr = divfunction.qubitnamesr
        qnd = divfunction.qubitnamesd
        qnc = divfunction.qubitnamec
        qnc = qnc[0]
        qn = qnq + qnr + qnd + [qnc]

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


if __name__ == "__main__":

    pass
