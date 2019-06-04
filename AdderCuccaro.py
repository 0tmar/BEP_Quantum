from cQASM import *


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


def cMAJ(qna, qnb, qnc, qnctrl):
    gates = []
    gates += [Qgate("toffoli", qnctrl, qna, qnb)]
    gates += [Qgate("cx", qna, qnc)]
    gates += [Qgate("toffoli", qnc, qnb, qna)]
    return gates


def cUMA(qna, qnb, qnc, qnctrl):
    gates = []
    gates += [Qgate("toffoli", qnc, qnb, qna)]
    gates += [Qgate("cx", qna, qnc)]
    gates += [Qgate("toffoli", qnctrl, qnc, qnb)]
    return gates


def COMPLEMENT_gates(n=None, qubitnames=None):

        qn = buildnames(n, qubitnames)

        gates = []

        for i in range(n):
            gates += [Qgate('x', qn[i])]

        return gates


def COMPLEMENT_gates_controlled(n=None, qubitnames=None, qubitnamectrl=None):

        qn = buildnames(n, qubitnames)
        qnctrl = buildnames(1, qubitnamectrl)

        gates = []

        for i in range(n):
            gates += [Qgate('cx', qnctrl[0], qn[i])]

        return gates


class COMPLEMENT(Qsubroutine):

    def __init__(self, n=None, qubitnames=None):

        qn = buildnames(n, qubitnames)
        self.qubitnames = qn

        gates = []

        for i in range(n):
            gates += [Qgate('x', qn[i])]

        super().__init__(name="complement", gates=gates)


class ADD(Qsubroutine):

    def __init__(self, n=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamez=None, do_overflow=True):

        if type(do_overflow) is not type(True):
            raise TypeError("'do_overflow' must be of type Boolean, not '{}'".format(type(do_overflow)))

        qna = buildnames(n, qubitnamesa, "a")
        qnb = buildnames(n, qubitnamesb, "b")
        if qubitnamec is None:
            qnc = "c"
        else:
            qnc = qubitnamec
        if do_overflow:
            if qubitnamez is None:
                qnz = "z"
            else:
                qnz = qubitnamez
        else:
            qnz = None

        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        self.qubitnamec = qnc
        if do_overflow:
            self.qubitnamez = qnz

        gates = []

        gates += MAJ(qna[0], qnb[0], qnc)
        for i in range(1, n):
            gates += MAJ(qna[i], qnb[i], qna[i-1])
        if do_overflow:
            gates += [Qgate("cx", qna[-1], qnz)]
        for i in range(n-1, 0, -1):
            gates += UMA(qna[i], qnb[i], qna[i-1])
        gates += UMA(qna[0], qnb[0], qnc)

        super().__init__(name="add", gates=gates)


class SUB(Qsubroutine):

    def __init__(self, n=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamez=None, do_overflow=True, subtype='a-b'):

        addsubroutine = ADD(n=n, qubitnamesa=qubitnamesa, qubitnamesb=qubitnamesb, qubitnamec=qubitnamec, qubitnamez=qubitnamez, do_overflow=do_overflow)
        addgates = addsubroutine.gates
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        qnc = addsubroutine.qubitnamec
        if do_overflow:
            qnz = addsubroutine.qubitnamez
        else:
            qnz = None
        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        self.qubitnamec = qnc
        if do_overflow:
            self.qubitnamez = qnz

        if subtype == "a-b":
            initcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qna)
            endcomplementgates = COMPLEMENT_gates(n=2*n, qubitnames=qna+qnb)
        elif subtype == "b-a":
            initcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qnb)
            endcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qnb)
        else:
            raise ValueError("Input of 'type' in cSUB function must either be 'a-b' of 'b-a', not '{}'".format(subtype))

        gates = initcomplementgates + addgates + endcomplementgates

        super().__init__(name="sub", gates=gates)


class cADD(Qsubroutine):

    def __init__(self, n=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamez=None, do_overflow=True, qubitnamectrl=None):

        if type(do_overflow) is not type(True):
            raise TypeError("'do_overflow' must be of type Boolean, not '{}'".format(type(do_overflow)))

        qna = buildnames(n, qubitnamesa, "a")
        qnb = buildnames(n, qubitnamesb, "b")
        if qubitnamec is None:
            qnc = "c"
        else:
            qnc = qubitnamec
        if do_overflow:
            if qubitnamez is None:
                qnz = "z"
            else:
                qnz = qubitnamez
        else:
            qnz = None
        if qubitnamectrl is None:
            qnctrl = "ctrl"
        else:
            qnctrl = qubitnamectrl

        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        self.qubitnamec = qnc
        if do_overflow:
            self.qubitnamez = qnz
        self.qubitnamectrl = qnctrl

        gates = []

        gates += cMAJ(qna[0], qnb[0], qnc, qnctrl)
        for i in range(1, n):
            gates += cMAJ(qna[i], qnb[i], qna[i-1], qnctrl)
        if do_overflow:
            gates += [Qgate("toffoli", qnctrl, qna[-1], qnz)]
        for i in range(n-1, 0, -1):
            gates += cUMA(qna[i], qnb[i], qna[i-1], qnctrl)
        gates += cUMA(qna[0], qnb[0], qnc, qnctrl)

        super().__init__(name="cadd", gates=gates)


class cSUB(Qsubroutine):

    def __init__(self, n=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamez=None, qubitnamectrl=None, do_overflow=True, subtype="a-b"):

        caddsubroutine = cADD(n=n, qubitnamesa=qubitnamesa, qubitnamesb=qubitnamesb, qubitnamec=qubitnamec, qubitnamez=qubitnamez, qubitnamectrl=qubitnamectrl, do_overflow=do_overflow)
        caddgates = caddsubroutine.gates
        qna = caddsubroutine.qubitnamesa
        qnb = caddsubroutine.qubitnamesb
        qnc = caddsubroutine.qubitnamec
        if do_overflow:
            qnz = caddsubroutine.qubitnamez
        else:
            qnz = None
        qnctrl = caddsubroutine.qubitnamectrl
        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        self.qubitnamec = qnc
        if do_overflow:
            self.qubitnamez = qnz
        self.qubitnamectrl = qnctrl

        if subtype == "a-b":
            initcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qna)
            endcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qna) + COMPLEMENT_gates_controlled(n=n, qubitnames=qna, qubitnamectrl=qnctrl)
        elif subtype == "b-a":
            initcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qnb)
            endcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qnb)
        else:
            raise ValueError("Input of 'type' in cSUB function must either be 'a-b' of 'b-a', not '{}'".format(subtype))

        gates = initcomplementgates + caddgates + endcomplementgates

        super().__init__(name="csub", gates=gates)


class ADDcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0", do_overflow=True):
        name = "Cuccaro Quantum Adder"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n + 2
        addsubroutine = ADD(n=n, do_overflow=do_overflow)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        qnc = addsubroutine.qubitnamec
        if do_overflow:
            qnz = addsubroutine.qubitnamez
        else:
            qnz = "z"

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

    def __init__(self, inp_a="0", inp_b="0", subtype="a-b", do_overflow=True):
        name = "Cuccaro Quantum Subtractor"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n + 2
        subsubroutine = SUB(n=n, subtype=subtype, do_overflow=do_overflow)
        qna = subsubroutine.qubitnamesa
        qnb = subsubroutine.qubitnamesb
        qnc = subsubroutine.qubitnamec
        if do_overflow:
            qnz = subsubroutine.qubitnamez
        else:
            qnz = "z"

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

        subroutines = [initsubroutine, subsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


class cADDcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0", inp_ctrl="0", do_overflow=True):
        name = "Controlled Cuccaro Quantum Adder"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n + 3
        addsubroutine = cADD(n=n, do_overflow=do_overflow)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        qnc = addsubroutine.qubitnamec
        if do_overflow:
            qnz = addsubroutine.qubitnamez
        else:
            qnz = "z"
        qnctrl = addsubroutine.qubitnamectrl

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        initgates += [Qgate("map", "q0", qnctrl)]
        initgates += [Qgate("map", "q1", qnc)]
        for i in range(n):
            initgates += [Qgate("map", "q"+str(2*i+2), qnb[i])]
            initgates += [Qgate("map", "q"+str(2*i+3), qna[i])]
        initgates += [Qgate("map", "q" + str(2*n + 2), qnz)]
        if inp_ctrl == "1":
            initgates += [Qgate("x", qnctrl)]
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


class cSUBcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0", inp_ctrl="0", subtype="a-b", do_overflow=True):
        name = "Controlled Cuccaro Quantum Subtractor"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n + 3
        subsubroutine = cSUB(n=n, subtype=subtype, do_overflow=do_overflow)
        qna = subsubroutine.qubitnamesa
        qnb = subsubroutine.qubitnamesb
        qnc = subsubroutine.qubitnamec
        if do_overflow:
            qnz = subsubroutine.qubitnamez
        else:
            qnz = "z"
        qnctrl = subsubroutine.qubitnamectrl

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        initgates += [Qgate("map", "q0", qnctrl)]
        initgates += [Qgate("map", "q1", qnc)]
        for i in range(n):
            initgates += [Qgate("map", "q"+str(2*i+2), qnb[i])]
            initgates += [Qgate("map", "q"+str(2*i+3), qna[i])]
        initgates += [Qgate("map", "q" + str(2*n + 2), qnz)]
        if inp_ctrl == "1":
            initgates += [Qgate("x", qnctrl)]
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

        subroutines = [initsubroutine, subsubroutine, resultsubroutine]

        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


if __name__ == "__main__":

    # add_circ = ADDcircuit(inp_a="11001", inp_b="00110")
    sub_circ = SUBcircuit(inp_a="11001", inp_b="00110")
    # print(add_circ)
    print(sub_circ)
