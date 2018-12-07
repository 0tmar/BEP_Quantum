from cQASM import *
from QFT import *

class ADD(Qsubroutine):

    def __init__(self, na=1, nb=1, qubitnamesa=None, qubitnamesb=None):

        if nb > na:
            raise ValueError("na must be greater than or equal to nb")

        qna = buildnames(na, qubitnamesa, "a")
        qnb = buildnames(nb, qubitnamesb, "b")

        self.qubitnamesa = qna
        self.qubitnamesb = qnb

        gates = []

        for i in range(int(nb)):
            j = na - nb + i
            gates += [Qgate('cz', qnb[i], qna[j])]
            if j is not 0:
                gates += [Qgate('swap', qnb[i], qna[j])]
                for k in range(j, 0, -1):
                    gates += [Qgate('cr', qna[j], qna[k])]
                gates += [Qgate('swap', qnb[i], qna[j])]
            # if i is not nb - 1:
            #     gates += [Qgate()]

        super().__init__(name="add", gates=gates)

class SUB(Qsubroutine):

    def __init__(self, na=1, nb=1, qubitnamesa=None, qubitnamesb=None):

        if nb > na:
            raise ValueError("na must be greater than or equal to nb")

        qna = buildnames(na, qubitnamesa, "a")
        qnb = buildnames(nb, qubitnamesb, "b")

        self.qubitnamesa = qna
        self.qubitnamesb = qnb

        gates = []

        for i in range(int(nb)):
            j = na - nb + i
            gates += [Qgate('cz', qnb[i], qna[j])]
            if j is not 0:
                gates += [Qgate('swap', qnb[i], qna[j])]
                for k in reversed(range(0, j)):
                    gates += [Qgate('cr', qna[j], qna[k])]
                    for m in range(j-1, k, -1):
                        gates += [Qgate('swap', qna[j], qna[m])]
                        gates += [Qgate('cr', qna[m], qna[k])]
                        gates += [Qgate('swap', qna[j], qna[m])]
                    gates += [Qgate('cz', qna[j], qna[k])]
                gates += [Qgate('swap', qnb[i], qna[j])]
            # if i is not nb - 1:
            #     gates += [Qgate()]

        super().__init__(name="sub", gates=gates)

class ADDcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0"):
        name = "Quantum Adder"
        na = len(inp_a)
        nb = len(inp_b)
        qubits = na + nb
        addsubroutine = ADD(na=na, nb=nb)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        for i in range(nb):
            initgates += [Qgate("map", "q"+str(i), qnb[i])]
        for i in range(na):
            initgates += [Qgate("map", "q"+str(nb+i), qna[i])]
        for i in range(na):
            if inp_a[i] == "1":
                initgates += [Qgate("x", qna[i])]
        for i in range(nb):
            if inp_b[i] == "1":
                initgates += [Qgate("x", qnb[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        qftsubroutine = QFT(n=na, qubitnames=qna)
        iqftsubroutine = iQFT(n=na, qubitnames=qna)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, qftsubroutine, addsubroutine, iqftsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)

class SUBcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0"):
        name = "Quantum Subtractor"
        na = len(inp_a)
        nb = len(inp_b)
        qubits = na + nb
        subsubroutine = SUB(na=na, nb=nb)
        qna = subsubroutine.qubitnamesa
        qnb = subsubroutine.qubitnamesb

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        for i in range(nb):
            initgates += [Qgate("map", "q"+str(i), qnb[i])]
        for i in range(na):
            initgates += [Qgate("map", "q"+str(nb+i), qna[i])]
        for i in range(na):
            if inp_a[i] == "1":
                initgates += [Qgate("x", qna[i])]
        for i in range(nb):
            if inp_b[i] == "1":
                initgates += [Qgate("x", qnb[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        qftsubroutine = QFT(n=na, qubitnames=qna)
        iqftsubroutine = iQFT(n=na, qubitnames=qna)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, qftsubroutine, subsubroutine, iqftsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)

if __name__ == "__main__":
    #add = ADD(na=5, nb=5)
    #print(add)

    add_circ = ADDcircuit(inp_a="11001", inp_b="00110")
    print(add_circ)
