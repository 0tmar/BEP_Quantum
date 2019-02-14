from cQASM import *
from QFT import QFT, iQFT
from AdderQFT import cADD, cSUB


def MUL(na, nb, nc, qubitnamesa=None, qubitnamesb=None, qubitnamesc=None, qubitnamez=None):

    if not nc >= na+nb:
        raise ValueError("nc must be greater than or equal to na+nb")

    qna = buildnames(n=na, qubitnames=qubitnamesa, defaultname="a")
    qnb = buildnames(n=nb, qubitnames=qubitnamesb, defaultname="b")
    qnc = buildnames(n=nc, qubitnames=qubitnamesc, defaultname="c")
    qnz = buildnames(n=1, qubitnames=qubitnamez, defaultname="z")

    subroutines = []

    for i in range(na):
        qftsubroutine = QFT(n=nc - i, qubitnames=qnc[:nc - i])
        addsubroutine = cADD(na=nc - i, nb=nb, qubitnamesa=qnc[:nc - i], qubitnamesb=qnb, qubitnamec=qna[na - 1 - i],
                             qubitnamed=qnz)
        iqftsubroutine = iQFT(n=nc - i, qubitnames=qnc[:nc - i])
        subroutines += [qftsubroutine, addsubroutine, iqftsubroutine]

    return subroutines


def MULSUB(na, nb, nc, qubitnamesa=None, qubitnamesb=None, qubitnamesc=None, qubitnamez=None):

    if not nc >= na+nb:
        raise ValueError("nc must be greater than or equal to na+nb")

    qna = buildnames(n=na, qubitnames=qubitnamesa, defaultname="a")
    qnb = buildnames(n=nb, qubitnames=qubitnamesb, defaultname="b")
    qnc = buildnames(n=nc, qubitnames=qubitnamesc, defaultname="c")
    if qubitnamez is None:
        qnz = ["z"]
    else:
        qnz = qubitnamez

    subroutines = []

    for i in range(na):
        qftsubroutine = QFT(n=nc - i, qubitnames=qnc[:nc - i])
        addsubroutine = cSUB(na=nc - i, nb=nb, qubitnamesa=qnc[:nc - i], qubitnamesb=qnb, qubitnamec=qna[na - 1 - i],
                             qubitnamed=qnz)
        iqftsubroutine = iQFT(n=nc - i, qubitnames=qnc[:nc - i])
        subroutines += [qftsubroutine, addsubroutine, iqftsubroutine]

    return subroutines


class MULcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0"):

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        name = "Cuccaro Quantum Multiplier"
        na = len(inp_a)
        nb = len(inp_b)
        nc = na+nb
        if na > nb:
            raise ValueError("len(a) must be smaller than or equal to len(b)")
        qubits = na + nb + nc + 1
        qna = buildnames(n=na, qubitnames="a")
        qnb = buildnames(n=nb, qubitnames="b")
        qnc = buildnames(n=na+nb, qubitnames="c")
        qnz = ["z"]

        qn = qna + qnb + qnc + qnz

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(na):
            if inp_a[i] == "1":
                initgates += [Qgate("x", qna[i])]
        for i in range(nb):
            if inp_b[i] == "1":
                initgates += [Qgate("x", qnb[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        subroutines = []
        subroutines += [initsubroutine]

        qftsubroutines_lst = []
        addsubroutines_lst = []
        iqftsubroutines_lst = []
        for i in range(na):
            qftsubroutine = QFT(n=nc-i, qubitnames=qnc[:nc-i])
            addsubroutine = cADD(na=nc-i, nb=nb, qubitnamesa=qnc[:nc-i], qubitnamesb=qnb, qubitnamec=qna[na-1-i], qubitnamed=qnz)
            iqftsubroutine = iQFT(n=nc-i, qubitnames=qnc[:nc-i])
            qftsubroutines_lst += [qftsubroutine]
            addsubroutines_lst += [addsubroutine]
            iqftsubroutines_lst += [iqftsubroutine]
            subroutines += [qftsubroutine, addsubroutine, iqftsubroutine]

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines += [resultsubroutine]

        super().__init__(name=name, qubits=qubits, subroutines=subroutines)
