from cQASM import *
from QFT import *


class ADD(Qsubroutine):

    def __init__(self, na=1, nb=1, qubitnamesa=None, qubitnamesb=None):
        """a: first number, to which b will be added (size na);
           b: second number, which will be added to a (size nb<=na);
           (a, b) -> (a+b, b)"""
        # Note: the first qubit (e.g. q0) is the MSB

        if nb > na:
            raise ValueError("len(a) must be greater than or equal to len(b)")

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
                for k in range(j-1, -1, -1):
                    gates += [Qgate('cr', qna[j], qna[k])]
                gates += [Qgate('swap', qnb[i], qna[j])]
            # if i is not nb - 1:
            #     gates += [Qgate()]

        super().__init__(name="add", gates=gates)


class SUB(Qsubroutine):

    def __init__(self, na=1, nb=1, qubitnamesa=None, qubitnamesb=None):
        """a: first number, from which b will be subtracted (size na);
           b: second number, which will be subtracted from a (size nb<=na);
           (a, b) -> (a-b, b)"""

        if nb > na:
            raise ValueError("len(a) must be greater than or equal to len(b)")

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


class cADD(Qsubroutine):

    def __init__(self, na=1, nb=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamed=None):
        """a: first number, to which b will be added (size na);
           b: second number, which will be added to a (size nb<=na);
           c: control whether the addition is performed (size 1);
           d: ancilla qubit (size 1);
           (a, b, c, d) -> (a-b, b, c, d)"""

        if nb > na:
            raise ValueError("len(a) must be greater than or equal to len(b)")

        qna = buildnames(na, qubitnamesa, "a")
        qnb = buildnames(nb, qubitnamesb, "b")
        qnc = buildnames(1, qubitnamec, "c")
        qnd = buildnames(1, qubitnamed, "d")

        qnc = qnc[0]
        qnd = qnd[0]

        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        self.qubitnamec = qnc
        self.qubitnamed = qnd

        gates = []

        for i in range(int(nb)):
            gates += [Qgate('toffoli', qnc, qnb[i], qnd)]
            j = na - nb + i
            gates += [Qgate('cz', qnd, qna[j])]
            if j is not 0:
                gates += [Qgate('swap', qnd, qna[j])]
                for k in range(j-1, -1, -1):
                    gates += [Qgate('cr', qna[j], qna[k])]
                gates += [Qgate('swap', qnd, qna[j])]
            gates += [Qgate('toffoli', qnc, qnb[i], qnd)]
            # if i is not nb - 1:
            #     gates += [Qgate()]

        super().__init__(name="cadd", gates=gates)


class cSUB(Qsubroutine):

    def __init__(self, na=1, nb=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamed=None):
        """a: first number, from which b will be subtracted (size na);
           b: second number, which will be subtracted from a (size nb<=na);
           c: control whether the subtraction is performed (size 1);
           d: ancilla qubit (size 1);
           (a, b, c, d) -> (a-b, b, c, d)"""

        if nb > na:
            raise ValueError("len(a) must be greater than or equal to len(b)")

        qna = buildnames(na, qubitnamesa, "a")
        qnb = buildnames(nb, qubitnamesb, "b")
        qnc = buildnames(1, qubitnamec, "c")
        qnd = buildnames(1, qubitnamed, "d")

        qnc = qnc[0]
        qnd = qnd[0]

        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        self.qubitnamec = qnc
        self.qubitnamed = qnd

        gates = []

        for i in range(int(nb)):
            gates += [Qgate('toffoli', qnc, qnb[i], qnd)]
            j = na - nb + i
            gates += [Qgate('cz', qnd, qna[j])]
            if j is not 0:
                gates += [Qgate('swap', qnd, qna[j])]
                for k in reversed(range(0, j)):
                    gates += [Qgate('cr', qna[j], qna[k])]
                    for m in range(j-1, k, -1):
                        gates += [Qgate('swap', qna[j], qna[m])]
                        gates += [Qgate('cr', qna[m], qna[k])]
                        gates += [Qgate('swap', qna[j], qna[m])]
                    gates += [Qgate('cz', qna[j], qna[k])]
                gates += [Qgate('swap', qnd, qna[j])]
            gates += [Qgate('toffoli', qnc, qnb[i], qnd)]
            # if i is not nb - 1:
            #     gates += [Qgate()]

        super().__init__(name="csub", gates=gates)


class ADDcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0"):
        name = "QFT Quantum Adder"
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
        name = "QFT Quantum Subtractor"
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


class cADDcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0", inp_c="0"):

        if (not isinstance(inp_a, str)) or (not isinstance(inp_b, str)) or (not isinstance(inp_c, str)):
            raise TypeError("input must be of type string")
        if not len(inp_c) == 1:
            raise ValueError("the length of input c must be equal to one")

        name = "Controlled QFT Quantum Adder"
        na = len(inp_a)
        nb = len(inp_b)
        qubits = na + nb + 2
        addsubroutine = cADD(na=na, nb=nb)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        qnc = addsubroutine.qubitnamec
        qnd = addsubroutine.qubitnamed

        initgates = []
        initgates += [Qgate("map", "q" + str(0), qnc)]
        for i in range(nb):
            initgates += [Qgate("map", "q"+str(i+1), qnb[i])]
        initgates += [Qgate("map", "q" + str(nb+1), qnd)]
        for i in range(na):
            initgates += [Qgate("map", "q"+str(nb+i+2), qna[i])]
        if inp_c == "1":
            initgates += [Qgate("x", qnc)]
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


class cSUBcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0", inp_c="0"):

        if (not isinstance(inp_a, str)) or (not isinstance(inp_b, str)) or (not isinstance(inp_c, str)):
            raise TypeError("input must be of type string")
        if not len(inp_c) == 1:
            raise ValueError("the length of input c must be equal to one")

        name = "Controlled QFT Quantum Subtractor"
        na = len(inp_a)
        nb = len(inp_b)
        qubits = na + nb + 2
        subsubroutine = cSUB(na=na, nb=nb)
        qna = subsubroutine.qubitnamesa
        qnb = subsubroutine.qubitnamesb
        qnc = subsubroutine.qubitnamec
        qnd = subsubroutine.qubitnamed

        initgates = []
        initgates += [Qgate("map", "q" + str(0), qnc)]
        for i in range(nb):
            initgates += [Qgate("map", "q"+str(i+1), qnb[i])]
        initgates += [Qgate("map", "q" + str(nb+1), qnd)]
        for i in range(na):
            initgates += [Qgate("map", "q"+str(nb+i+2), qna[i])]
        if inp_c == "1":
            initgates += [Qgate("x", qnc)]
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
    # add = ADD(na=5, nb=5)
    # print(add)

    # add_circ = ADDcircuit(inp_a="11001", inp_b="00110")
    cadd_circ = cADDcircuit(inp_a="11001", inp_b="00110")
    print(cadd_circ)
