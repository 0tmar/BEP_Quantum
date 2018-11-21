from cQASM import *

class QFT(Qsubroutine):

    def __init__(self, n=1, qubitnames=None):

        if isinstance(qubitnames, list):
            if len(qubitnames) == n:
                qn = qubitnames
            else:
                raise ValueError("Incorrect number of qubit names")
        elif isinstance(qubitnames, str):
            qn = []
            for i in range(n):
                qn += [qubitnames + str(i)]
        else:
            qn = []
            string = "q"
            for i in range(n):
                qn += [string + str(i)]

        self.qubitnames = qn

        gates = []

        for i in range(n):
            gates += [Qgate('h',[qn[i]])]
            for j in range(i+1,n):
                gates += [Qgate('cr', [qn[i],qn[j]])]

        for i in range(int(n/2)):
            gates += [Qgate('swap',[qn[i],qn[n-1-i]])]

        super().__init__(name="qft", gates=gates)

class QFTcircuit(Qfunc):

    def __init__(self, n=1, input="0"):
        name = "Quantum Fourier Transform"
        qubits = n
        qftsubroutine = QFT(n=n)
        qn = qftsubroutine.qubitnames

        if len(input) is not n:
            raise ValueError("Error must be of length n")
        if not isinstance(input, str):
            raise TypeError("input must be of type string")

        initgates = []
        for i in range(n):
            if input[i] == "1":
                initgates += [Qgate("x", [qn[i]])]
        initgates += [Qgate("display", [])]

        initsubroutine = Qsubroutine(name="init", gates=initgates)

        resultgates = [Qgate("measure", []), Qgate("display", [])]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = [initsubroutine, qftsubroutine, resultsubroutine]
        super().__init__(name=name, qubits=qubits, subroutines=subroutines)

if __name__ == "__main__":
    qft = QFT(n=10)
    #print(qft)

    qft_circ = QFTcircuit(n=5, input="11011")
    print(qft_circ)