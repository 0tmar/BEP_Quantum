from cQASM import *


class COMPLEMENT(Qsubroutine):

    def __init__(self, n=None, qubitnames=None):

        qn = buildnames(n, qubitnames)
        self.qubitnames = qn

        gates = []

        for i in range(n):
            gates += [Qgate('x', qn[i])]

        super().__init__(name="complement", gates=gates)


class ADD(Qsubroutine):

    def __init__(self, n=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamez=None):

        def MAJ(qna, qnb, qnc):
            gates = []
            gates += [Qgate("cx", qna, qnb)]
            gates += [Qgate("cx", qna, qnc)]
            gates += [Qgate("toffoli", qnc, qnb, qna)]
            return gates

        def UMA(qna, qnb, qnc):
            gates = []
            gates += [Qgate("toffoli", qnc, qnb, qna)]
            gates += [Qgate("cx", qna, qnc)]
            gates += [Qgate("cx", qnc, qnb)]
            return gates

        qna = buildnames(n, qubitnamesa, "a")
        qnb = buildnames(n, qubitnamesb, "b")
        if qubitnamez == None:
            qnc = "c"
        else:
            qnc = qubitnamec
        if qubitnamez == None:
            qnz = "z"
        else:
            qnz = qubitnamez

        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        self.qubitnamec = qnc
        self.qubitnamez = qnz

        gates = []

        gates += MAJ(qna[0], qnb[0], qnc)
        for i in range(1, n):
            gates += MAJ(qna[i], qnb[i], qna[i-1])
        gates += [Qgate("cx", qna[-1], qnz)]
        for i in range(n-1, 0, -1):
            gates += UMA(qna[i], qnb[i], qna[i-1])
        gates += UMA(qna[0], qnb[0], qnc)

        super().__init__(name="add", gates=gates)


class ADDcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0"):
        name = "Quantum Adder"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n + 2
        addsubroutine = ADD(n=n)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        qnc = addsubroutine.qubitnamec
        qnz = addsubroutine.qubitnamez

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        initgates += [Qgate("map", "q0", qnc)]
        for i in range(n):
            initgates += [Qgate("map", "q"+str(2*i+1), qnb[i])]
            initgates += [Qgate("map", "q"+str(2*i+2), qna[i])]
        initgates += [Qgate("map", "q" + str(2*n + 1), qnz)]
        for i in range(n):
            if inp_a[-i-1] == "1":
                initgates += [Qgate("x", qna[i])]
        for i in range(n):
            if inp_b[-i-1] == "1":
                initgates += [Qgate("x", qnb[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, addsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


class SUBcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0"):
        name = "Quantum Subtractor"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n + 2
        addsubroutine = ADD(n=n)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        qnc = addsubroutine.qubitnamec
        qnz = addsubroutine.qubitnamez

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        initgates += [Qgate("map", "q0", qnc)]
        for i in range(n):
            initgates += [Qgate("map", "q"+str(2*i+1), qnb[i])]
            initgates += [Qgate("map", "q"+str(2*i+2), qna[i])]
        initgates += [Qgate("map", "q" + str(2*n + 1), qnz)]
        for i in range(n):
            if inp_a[-i-1] == "1":
                initgates += [Qgate("x", qna[i])]
        for i in range(n):
            if inp_b[-i-1] == "1":
                initgates += [Qgate("x", qnb[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        initcomplementsubroutine = COMPLEMENT(n=n, qubitnames=qna)
        endcomplementsubroutine = COMPLEMENT(n=2*n, qubitnames=qna+qnb)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, initcomplementsubroutine, addsubroutine, endcomplementsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


if __name__ == "__main__":

    # add_circ = ADDcircuit(inp_a="11001", inp_b="00110")
    sub_circ = SUBcircuit(inp_a="11001", inp_b="00110")
    # print(add_circ)
    print(sub_circ)