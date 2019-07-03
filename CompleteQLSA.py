import math
from cQASM import *
from QFT import *
from Cao2012Experiment import *
from DivisionThapliyal import *
from AncillaRotation import *


class CompleteQlsaWithCaoMatrix(Qfunction):

    def __init__(self, n_eig_inv, k_anc_rot, r_anc_rot, input_eigenvector=None, input_n_state=None):

        n_vec_b = 2   # Number of qubits used for saving vector b
        n_lambda = 4  # Number of qubits used to store the eigenvalues

        if n_eig_inv < n_lambda:
            raise ValueError("Value for n_eig_inv must be equal to or greater than n_lambda = 4.")

        # 2 qubits for storing |b>, 4 qubits for storing the eigenvalues, (4+ n_eig_inv + 2) for inverting them,
        # and (1 + max(0, min(n_eig_inv+1, 1+2*k_anc_rot) - 2)) qubits for the ancilla rotation.
        qubits = n_vec_b + n_lambda + (n_lambda + n_eig_inv + 3) + (1 + max(0, min(n_eig_inv, 1 + 2 * k_anc_rot) - 2))

        qna = buildnames(n=n_vec_b, qubitnames="a")                # qubits storing vector |b>
        qnb = buildnames(n=n_lambda+1, qubitnames="b")             # qubits storing the eigenvalues
        qnc = buildnames(n=n_eig_inv+1, qubitnames="c")            # qubits storing the inverted eigenvalues
        qnd = buildnames(n=n_lambda+1, qubitnames="d")             # ancilla qubits for inverting the eigenvalues
        qne = buildnames(n=max(0, min(n_eig_inv, 1 + 2 * k_anc_rot) - 2), qubitnames="e")  # ancilla qubits for the ancilla rotation
        qnf = buildnames(n=1, qubitnames="f")                      # ancilla qubit which is being rotated

        qn = qnf + qne + qnd + qnc + qnb + qna

        initgates = []
        for i in range(qubits):
            initgates += [Qgate("map", "q"+str(i), qn[i])]
        for i in range(n_lambda):
            initgates += [Qgate("h", qnb[i])]
        if input_n_state is None:
            for i in range(n_vec_b):
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
        initgates += [Qgate("x", qnd[-1])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        if n_eig_inv == n_lambda:
            ancrotgates = list(reversed(qnc[1:]))
        elif n_eig_inv == n_lambda + 1:
            ancrotgates = list(reversed(qnc))
        else:
            ancrotgates = list(reversed(qnd[(-n_eig_inv + n_lambda + 1):] + qnc))

        expAsubroutines = []
        for i in range(n_lambda):
            expAsubroutines += [expA(
                qubitnamea=qnb[i],
                qubitnameb=qna[0],
                qubitnamec=qna[1],
                sign=1,
                n=i)]

        iqftsubroutine = iQFT(
            n=n_lambda,
            qubitnames=qnb[:-1])

        reversesubroutine = REVERSE(
            n=n_lambda,
            qubitnames=qnb[:-1])

        divsubroutine = DIV(
            n=n_eig_inv+1,
            m=n_lambda+1,
            qubitnamesn=qnd,
            qubitnameso=qnc,
            qubitnamesd=qnb,
            sign=1)

        ancrotsubroutine = AncillaRotation(
            n=n_eig_inv,  # In general this should be n_eig_inv+1, but it is not necessary due to the guaranteed zero
            c=math.pi/(2**(r_anc_rot-1)),
            k=k_anc_rot,
            qubitnamesc=ancrotgates,
            qubitnamer=qnf,
            qubitnamesa=qne)

        undivsubroutine = DIV(
            n=n_eig_inv+1,
            m=n_lambda+1,
            qubitnamesn=qnd,
            qubitnameso=qnc,
            qubitnamesd=qnb,
            sign=-1)

        qftsubroutine = QFT(
            n=n_lambda,
            qubitnames=qnb[:-1])

        displaysubroutine = Qsubroutine(name="display", gates=[Qgate("display")])

        unexpAsubroutines = []
        for i in reversed(range(n_lambda)):
            unexpAsubroutines += [expA(
                qubitnamea=qnb[i],
                qubitnameb=qna[0],
                qubitnamec=qna[1],
                sign=-1,
                n=i)]

        uninitgates = []
        for i in range(n_lambda):
            uninitgates += [Qgate("h", qnb[i])]
        uninitgates += [Qgate("x", qnd[-1])]
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
