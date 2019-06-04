from cQASM import *


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


class ADD(Qsubroutine):

    def __init__(self, n=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamez=None, qubitnamectrl=None, do_overflow=True, do_ctrl=True):

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
        if do_ctrl:
            if qubitnamectrl is None:
                qnctrl = "ctrl"
            else:
                qnctrl = qubitnamectrl
        else:
            qnctrl = None

        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        if do_overflow and do_ctrl:
            self.qubitnamec = qnc
        if do_overflow:
            self.qubitnamez = qnz
        if do_ctrl:
            self.qubitnamectrl = qnctrl

        gates = []

        for i in range(n-1):
            gates += [Qgate("cx", qna[i+1], qnb[i+1])]
        # if do_overflow:
        #     if do_ctrl:
        #         gates += [Qgate("toffoli", qna[-1], qnb[-1], qnc)]
        #         gates += [Qgate("toffoli", qnctrl, qnc, qnz)]
        #         gates += [Qgate("toffoli", qna[-1], qnb[-1], qnc)]
        #     else:
        #         gates += [Qgate("toffoli", qna[-1], qnb[-1], qnz)]
        for i in range(n-2):
            gates += [Qgate("cx", qna[-i-2], qna[-i-1])]
        for i in range(n-1):
            gates += [Qgate("toffoli", qna[i], qnb[i], qna[i+1])]
        if do_overflow:
            if do_ctrl:
                gates += [Qgate("toffoli", qna[-1], qnb[-1], qnc)]
                gates += [Qgate("display")]
                gates += [Qgate("toffoli", qnctrl, qnc, qnz)]
                gates += [Qgate("display")]
                gates += [Qgate("toffoli", qna[-1], qnb[-1], qnc)]
            else:
                gates += [Qgate("toffoli", qna[-1], qnb[-1], qnz)]
        for i in range(n-1):
            if do_ctrl:
                gates += [Qgate("toffoli", qnctrl, qna[-i-1], qnb[-i-1])]
            else:
                gates += [Qgate("cx", qna[-i-1], qnb[-i-1])]
            gates += [Qgate("toffoli", qna[-i-2], qnb[-i-2], qna[-i-1])]
        if do_ctrl:
            gates += [Qgate("toffoli", qnctrl, qna[0], qnb[0])]
        else:
            gates += [Qgate("cx", qna[0], qnb[0])]
        for i in range(n-2):
            gates += [Qgate("cx", qna[i+1], qna[i+2])]
        for i in range(n-1):
            gates += [Qgate("cx", qna[-i-1], qnb[-i-1])]

        super().__init__(name="add", gates=gates)


class SUB(Qsubroutine):

    def __init__(self, n=1, qubitnamesa=None, qubitnamesb=None, qubitnamec=None, qubitnamez=None, qubitnamectrl=None, do_overflow=True, do_ctrl=True, subtype="a-b"):

        addsubroutine = ADD(n=n, qubitnamesa=qubitnamesa, qubitnamesb=qubitnamesb, qubitnamec=qubitnamec, qubitnamez=qubitnamez, qubitnamectrl=qubitnamectrl, do_overflow=do_overflow, do_ctrl=do_ctrl)
        addgates = addsubroutine.gates
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        if do_overflow and do_ctrl:
            qnc = addsubroutine.qubitnamec
        else:
            qnc = None
        if do_overflow:
            qnz = addsubroutine.qubitnamez
        else:
            qnz = None
        if do_ctrl:
            qnctrl = addsubroutine.qubitnamectrl
        else:
            qnctrl = None
        self.qubitnamesa = qna
        self.qubitnamesb = qnb
        if do_overflow and do_ctrl:
            self.qubitnamec = qnc
        if do_overflow:
            self.qubitnamez = qnz
        if do_ctrl:
            self.qubitnamectrl = qnctrl

        if subtype == "a-b":
            initcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qna)
            if do_ctrl:
                endcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qna)
                if do_overflow:
                    endcomplementgates += COMPLEMENT_gates_controlled(n=n+1, qubitnames=qnb+[qnz], qubitnamectrl=qnctrl)
                else:
                    endcomplementgates += COMPLEMENT_gates_controlled(n=n, qubitnames=qnb, qubitnamectrl=qnctrl)
            else:
                endcomplementgates = COMPLEMENT_gates(n=2*n, qubitnames=qna+qnb)
        elif subtype == "b-a":
            if do_overflow:
                initcomplementgates = COMPLEMENT_gates(n=n+1, qubitnames=qnb+[qnz])
                endcomplementgates = COMPLEMENT_gates(n=n+1, qubitnames=qnb+[qnz])
            else:
                initcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qnb)
                endcomplementgates = COMPLEMENT_gates(n=n, qubitnames=qnb)
        else:
            raise ValueError("Input of 'type' in cSUB function must either be 'a-b' of 'b-a', not '{}'".format(subtype))

        gates = initcomplementgates + addgates + endcomplementgates

        super().__init__(name="sub", gates=gates)


class ADDcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0", inp_ctrl="0", do_overflow=True, do_ctrl=True):
        name = "Munoz-Coreas Quantum Adder"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n
        if do_overflow:
            qubits += 1
        if do_ctrl:
            qubits += 1
        if do_overflow and do_ctrl:
            qubits += 1
        addsubroutine = ADD(n=n, do_overflow=do_overflow, do_ctrl=do_ctrl)
        qna = addsubroutine.qubitnamesa
        qnb = addsubroutine.qubitnamesb
        if do_overflow and do_ctrl:
            qnc = addsubroutine.qubitnamec
        else:
            qnc = None
        if do_overflow:
            qnz = addsubroutine.qubitnamez
        else:
            qnz = "z"
        if do_ctrl:
            qnctrl = addsubroutine.qubitnamectrl
        else:
            qnctrl = None

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        j = 0
        if do_ctrl:
            initgates += [Qgate("map", "q0", qnctrl)]
            j += 1
        for i in range(n):
            initgates += [Qgate("map", "q"+str(2*i+j), qnb[i])]
            initgates += [Qgate("map", "q"+str(2*i+j+1), qna[i])]
        if do_overflow:
            initgates += [Qgate("map", "q" + str(2*n+j), qnz)]
        if do_overflow and do_ctrl:
            initgates += [Qgate("map", "q" + str(2*n+j), qnc)]
            j += 1
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


class SUBcircuit(Qfunction):

    def __init__(self, inp_a="0", inp_b="0", inp_ctrl="0", do_overflow=True, do_ctrl=True):
        name = "Munoz-Coreas Quantum Subtractor"
        na = len(inp_a)
        nb = len(inp_b)
        n = max(na, nb)
        inp_a = (n-na)*"0" + inp_a
        inp_b = (n-nb)*"0" + inp_b
        qubits = 2*n
        if do_overflow:
            qubits += 1
        if do_ctrl:
            qubits += 1
        if do_overflow and do_ctrl:
            qubits += 1
        subsubroutine = SUB(n=n, do_overflow=do_overflow, do_ctrl=do_ctrl)
        qna = subsubroutine.qubitnamesa
        qnb = subsubroutine.qubitnamesb
        if do_overflow and do_ctrl:
            qnc = subsubroutine.qubitnamec
        else:
            qnc = None
        if do_overflow:
            qnz = subsubroutine.qubitnamez
        else:
            qnz = "z"
        if do_ctrl:
            qnctrl = subsubroutine.qubitnamectrl
        else:
            qnctrl = None

        if not isinstance(inp_a, str):
            raise TypeError("input must be of type string")
        if not isinstance(inp_b, str):
            raise TypeError("input must be of type string")

        initgates = []
        j = 0
        if do_ctrl:
            initgates += [Qgate("map", "q0", qnctrl)]
            j += 1
        for i in range(n):
            initgates += [Qgate("map", "q"+str(2*i+j), qnb[i])]
            initgates += [Qgate("map", "q"+str(2*i+j+1), qna[i])]
        if do_overflow and do_ctrl:
            initgates += [Qgate("map", "q" + str(2*n+j), qnc)]
            j += 1
        if do_overflow:
            initgates += [Qgate("map", "q" + str(2*n+j), qnz)]
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
