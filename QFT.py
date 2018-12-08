from cQASM import *


class REVERSE(Qsubroutine):

    def __init__(self, n=1, qubitnames=None):

        qn = buildnames(n, qubitnames)
        self.qubitnames = qn

        gates = []

        for i in range(int(n/2)):
            gates += [Qgate('swap', qn[i], qn[n-1-i])]

        super().__init__(name="reverse", gates=gates)


class QFT(Qsubroutine):

    def __init__(self, n=1, qubitnames=None):

        qn = buildnames(n, qubitnames)
        self.qubitnames = qn

        gates = []

        for i in range(n):
            gates += [Qgate('h', qn[i])]
            for j in range(i+1,n):
                gates += [Qgate('cr', qn[j], qn[i])]

        super().__init__(name="qft", gates=gates)


class iQFT(Qsubroutine):

    def __init__(self, n=1, qubitnames=None):

        qn = buildnames(n, qubitnames)
        self.qubitnames = qn

        gates = []

        for i in reversed(range(n)):
            #gates += [Qgate('')]
            for j in reversed(range(i+1, n)):
                #gates += [Qgate('')]
                gates += [Qgate('cr', qn[j], qn[i])]
                for k in range(1, j-i):
                    gates += [Qgate('swap', qn[j], qn[j-k])]
                    gates += [Qgate('cr', qn[j-k], qn[i])]
                    gates += [Qgate('swap', qn[j], qn[j-k])]
                gates += [Qgate('cz', qn[j], qn[i])]
            gates += [Qgate('h', qn[i])]

        super().__init__(name="iqft", gates=gates)


class QFTcircuit(Qfunction):

    def __init__(self, inp="0"):
        name = "Quantum Fourier Transform"
        n = len(inp)
        qubits = n
        qftsubroutine = QFT(n=n)
        qn = qftsubroutine.qubitnames

        if not isinstance(inp, str):
            raise TypeError("input must be of type string")

        initgates = []
        for i in range(n):
            if inp[i] == "1":
                initgates += [Qgate("x", qn[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)
        reversesubroutine = REVERSE(n=n, qubitnames=qn)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, qftsubroutine, reversesubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


class QFT_iQFTcircuit(Qfunction):

    def __init__(self, inp="0"):
        name = "Back and forth Quantum Fourier Transform"
        n = len(inp)
        qubits = n
        qftsubroutine = QFT(n=n)
        qn = qftsubroutine.qubitnames

        if not isinstance(inp, str):
            raise TypeError("input must be of type string")

        initgates = []
        for i in range(n):
            if inp[i] == "1":
                initgates += [Qgate("x", qn[i])]
        initgates += [Qgate("display")]

        initsubroutine = Qsubroutine(name="init", gates=initgates)
        iqftsubroutine = iQFT(n=n, qubitnames=qn)

        resultgates = [Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, qftsubroutine, iqftsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)


if __name__ == "__main__":
    qft = QFT(n=10)
    # print(qft)

    qft_circ = QFT_iQFTcircuit(inp="11001")
    print(qft_circ)
