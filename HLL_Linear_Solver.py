from cQASM import *
from cQASM_Compound_Gates import *
from QFT import QFT, iQFT, REVERSE
from math import pi


def cRy(qubitnamea, qubitnameb, n, r):
    """"outputs controlled-Ry(2^(n)*pi/2^(r-1)) on qubit b, with qubit a as control"""

    qna = qubitnamea
    qnb = qubitnameb

    theta = pi*2**(n-r+1)

    gates = []
    gates += [Qgate('#', " performs Ry(2^({})*pi/2^({}-1)) = Ry(pi/{}) on {}, controlled by {}"
                         .format(n, r, 2**(r-n-1), qnb, qna))]
    gates += [Qgate('ry', qnb, theta/4)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('ry', qnb, -theta/2)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('ry', qnb, theta/4)]

    cRysubroutine = Qsubroutine(name="cRy", gates=gates)

    return cRysubroutine


def eigenvalue_inversion(qnl, qnm, qnc):
    """outputs the circuit part designed to invert the eigenvalues"""

    def crzz(qnl, qnm):
        gates = []
        gates += [Qgate('#', ' Performs cRzz(2*pi/2^(m-i)) on L, controlled by the i-th qubit of M, for all i in M')]
        for i in range(len(qnm)):
            # applies a cRzz(2*pi/2^(m-i)) to the last L qubit, controlled by the i-th M qubit
            # [= Rz(2*pi/2^(m-i)) on the i-th M qubit, up to global phase]
            gates += [Qgate('rz', qnm[i], 2*pi/(2**(len(qnm)-i)))]
        subroutine = Qsubroutine(name='cRzz', gates=gates)
        return subroutine

    def exp_iH0t2m(qnl, qnm, qnc, t):
        """outputs the gates to perform cc-exp(i*H0*t/2^m), which is ccRz(t/2^i) on the i-th qubit of M (for all i)"""
        # Watch out: qnl and qnc should be a specific qubit, while qnm should be all qubits of M!

        gates = []
        gates += [Qgate()]
        gates += [Qgate('#', '## Performs exp(i*H0*t/2^m) with t = {}'.format(t))]
        for i in range(len(qnm)):
            gates += ccRz(qna=qnl, qnb=qnc, qnc=qnm[i], theta=t/(2**(i+1)))
        return gates

    def G(qnl, qnm, qnc, l, kl, t):
        """outputs the gates to perform cc-G(l - kl), which is
        exp_iH0t2m(qnl[kl], qnm, qnc[kc], t=t0/2**(u+v-m)) for all kc in C (with u=c-kc and v=l-kl)"""
        # Watch out: qnl should be a specific qubit, while qnm and qnc should be all qubits of M and C!

        v = l - kl
        m = len(qnm)
        c = len(qnc)

        gates = []
        gates += [Qgate()]
        gates += [Qgate('#', '#### Performs G(l - kl) with l = {} and kl = {}'.format(l, kl))]
        for kc in range(len(qnc)):
            u = c - kc
            gates += exp_iH0t2m(qnl=qnl, qnm=qnm, qnc=qnc[kc], t=t/(2**(u + v - m)))
        return gates

    def exp_iH0t0(qnl, qnm, qnc):
        """outputs the gates to perform exp(i*H0*t0), which is (qnl[kl], qnm, qnc, l, kl, t0) for all kl in L
        with t0=2*pi"""

        t0=2*pi
        l = len(qnl)

        gates = []
        gates += [Qgate()]
        gates += [Qgate('#', '###### Performs exp(i*H0*t0) with t0 = 2*pi')]
        for kl in range(len(qnc)):
            gates += G(qnl[kl], qnm, qnc, l, kl, t0)
        return gates

    crzzsubroutine = crzz(qnl=qnl, qnm=qnm)

    expiH0t0gates = exp_iH0t0(qnl=qnl, qnm=qnm, qnc=qnc)
    expiH0t0subroutine = Qsubroutine(name='exp_iH0t0', gates=expiH0t0gates)

    subroutines = [crzzsubroutine, expiH0t0subroutine]

    return subroutines


class EigenvalueInversion_Circuit(Qfunction):

    def __init__(self, n=4, x=2):
        # n: size of C
        # x: which qubit of C should be flipped

        qnl = buildnames(n, qubitnames="l")
        qnm = buildnames(n + 4, qubitnames="m")
        qnc = buildnames(n, qubitnames="c")
        qn = qnl + qnm + qnc
        qndef = buildnames(n=len(qn))

        initgates = []
        for i in range(len(qn)):
            initgates += [Qgate('map', qndef[i], qn[i])]
        for i in range(len(qnl)):
            initgates += [Qgate('h', qnl[i])]
        for i in range(len(qnm)):
            initgates += [Qgate('h', qnm[i])]
        initgates += [Qgate('x', qnc[x])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        eigenvalueinversionsubroutine = eigenvalue_inversion(qnl=qnl, qnm=qnm, qnc=qnc)

        measuregates = [Qgate('display'), Qgate('measure')]
        measuresubroutine = Qsubroutine(name='measure', gates=measuregates)

        subroutines = [initsubroutine] + eigenvalueinversionsubroutine + [measuresubroutine]

        super().__init__(name='Eigenvalue inversion test', qubits=len(qn), subroutines=subroutines)


if __name__ == "__main__":
    circ = EigenvalueInversion_Circuit(n=4, x=2)
    print(circ)
