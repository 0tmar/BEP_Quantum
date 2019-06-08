import math
from cQASM import *
from QFT import *
from Cao2012_Experiment import *
from DivisionThapliyal import *
from Ancilla_Rotation import *


class CompleteQlsaWithCaoMatrix(Qfunction):

    def __init__(self, n_eig_inv, m_anc_rot, r_anc_rot, input_eigenvector=None, input_n_state=None):

        # 2 qubits for storing |b>, 4 qubits for storing the eigenvalues, (4+ n_eig_inv + 2) for inverting them,
        # and (1 + max(0, min(n_eig_inv+1, 1+2*m_anc_rot) - 2)) qubits for the ancilla rotation.
        qubits = 2 + 4 + (4 + n_eig_inv + 3) + (1 + max(0, min(n_eig_inv+1, 1+2*m_anc_rot) - 2))

        qna = buildnames(n=2, qubitnames="a")                      # qubits storing vector |b>
        qnb = buildnames(n=4+1, qubitnames="b")                    # qubits storing the eigenvalues
        qnc = buildnames(n=n_eig_inv+1, qubitnames="c")            # qubits storing the inverted eigenvalues
        qnd = buildnames(n=4+1, qubitnames="d")                    # ancilla qubits for inverting the eigenvalues
        qne = buildnames(n=max(0, min(n_eig_inv+1, 1+2*m_anc_rot) - 2), qubitnames="e")  # ancilla qubits for the ancilla rotation
        qnf = buildnames(n=1, qubitnames="f")                      # ancilla qubit which is being rotated

        qn = qna + qnb + qnc + qnd + qne + qnf

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(len(qnb)-1):
            initgates += [Qgate("h", qnb[i])]
        if input_n_state is None:
            for i in range(len(qna)):
                initgates += [Qgate("h", qna[i])]
        if input_eigenvector is not None:
            if input_eigenvector == 0 or input_eigenvector == 1:
                initgates += [Qgate("x", qna[0])]
            else:
                initgates += [Qgate("#", "x {}".format(qna[0]))]
            if input_eigenvector == 0 or input_eigenvector == 2:
                initgates += [Qgate("x", qna[1])]
            else:
                initgates += [Qgate("#", "x {}".format(qna[1]))]
            initgates += [Qgate("cz", qna[0], qna[1])]
            if input_eigenvector == 0 or input_eigenvector == 2:
                initgates += [Qgate("x", qna[1])]
            else:
                initgates += [Qgate("#", "x {}".format(qna[1]))]
            if input_eigenvector == 0 or input_eigenvector == 1:
                initgates += [Qgate("x", qna[0])]
            else:
                initgates += [Qgate("#", "x {}".format(qna[0]))]
        if input_n_state is not None:
            if input_n_state == 2 or input_n_state == 3:
                initgates += [Qgate("x", qna[0])]
            if input_n_state == 1 or input_n_state == 3:
                initgates += [Qgate("x", qna[1])]
        initgates += [Qgate("x", qnc[-1])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        expAsubroutines = []
        for i in range(len(qnb)-1):
            expAsubroutines += [expA(qubitnamea=qnb[i], qubitnameb=qna[0], qubitnamec=qna[1], sign=1, n=i)]

        iqftsubroutine = iQFT(n=len(qnb)-1, qubitnames=qnb[:-1])
        reversesubroutine = REVERSE(n=len(qnb)-1, qubitnames=qnb[:-1])
        divsubroutine = DIV(n=n_eig_inv+1, m=len(qnb), qubitnamesq=qnc, qubitnamesr=qnd, qubitnamesd=qnb, sign=1)
        ancrotsubroutine = AncillaRotation(n=n_eig_inv+1, c=math.pi/(2**(r_anc_rot-1)), m=m_anc_rot, qubitnamesc=list(reversed(qnc[(-n_eig_inv-1+len(qnd)):]+qnd)), qubitnamer=qnf, qubitnamesa=qne)
        undivsubroutine = DIV(n=n_eig_inv+1, m=len(qnb), qubitnamesq=qnc, qubitnamesr=qnd, qubitnamesd=qnb, sign=-1)
        qftsubroutine = QFT(n=len(qnb)-1, qubitnames=qnb[:-1])

        displaysubroutine = Qsubroutine(name="display", gates=[Qgate("display")])

        unexpAsubroutines = []
        for i in reversed(range(len(qnb)-1)):
            unexpAsubroutines += [expA(qubitnamea=qnb[i], qubitnameb=qna[0], qubitnamec=qna[1], sign=-1, n=i)]

        uninitgates = []
        for i in range(len(qnb)-1):
            uninitgates += [Qgate("h", qnb[i])]
        uninitgates += [Qgate("x", qnc[-1])]
        uninitsubroutine = Qsubroutine(name="uninit", gates=uninitgates)

        resultgates = [Qgate("display"), Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = []
        subroutines += [initsubroutine]
        subroutines += expAsubroutines
        subroutines += [iqftsubroutine]
        subroutines += [reversesubroutine]
        subroutines += [displaysubroutine]
        subroutines += [divsubroutine]
        subroutines += [displaysubroutine]
        subroutines += [ancrotsubroutine]
        subroutines += [undivsubroutine]
        subroutines += [reversesubroutine]
        subroutines += [qftsubroutine]
        subroutines += unexpAsubroutines
        subroutines += [uninitsubroutine]
        subroutines += [resultsubroutine]

        super().__init__(name="Complete QLSA for the Cao matrix", qubits=qubits, subroutines=subroutines)
