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
            if i is not na-1:
                gates += [Qgate('swap', qnb[i], qna[na-j-2])]
            for k in range(na-j-1,na):
                gates += [Qgate('cr', qnb[i], qna[k])]
            if i is not na-1:
                gates += [Qgate('swap', qnb[i], qna[na-j-2])]
            #if i is not nb - 1:
            #    gates += [Qgate()]

        super().__init__(name="add", gates=gates)

class ADDcircuit(Qfunction):

    def __init__(self, inputa="0", inputb="0"):
        name = "Quantum Adder"
        na = len(inputa)
        nb = len(inputb)
        qubits = na + nb
        addsubroutine = ADD(na=na, nb=nb)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb

        if not isinstance(inputa, str):
            raise TypeError("input must be of type string")
        if not isinstance(inputb, str):
            raise TypeError("input must be of type string")

        initgates = []
        for i in range(nb):
            initgates += [Qgate("map", "q"+str(i), qnb[i])]
        for i in range(na):
            initgates += [Qgate("map", "q"+str(nb+i), qna[i])]
        for i in range(na):
            if inputa[i] == "1":
                initgates += [Qgate("x", qna[i])]
        for i in range(nb):
            if inputb[i] == "1":
                initgates += [Qgate("x", qnb[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        qftsubroutine = QFT(n=na, qubitnames=qna)
        iqftsubroutine = iQFT(n=na, qubitnames=qna)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, qftsubroutine, addsubroutine, iqftsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)

if __name__ == "__main__":
    #add = ADD(na=5, nb=5)
    #print(add)

    add_circ = ADDcircuit(inputa="11011", inputb="00100")
    print(add_circ)
