from cQASM import *
from cQASM_Compound_Gates import *
from QFT import QFT, iQFT, REVERSE
from math import pi


# def cRy(qubitnamea, qubitnameb, n, r):
#     """"outputs controlled-Ry(2^(n)*pi/2^(r-1)) on qubit b, with qubit a as control"""
#
#     qna = qubitnamea
#     qnb = qubitnameb
#
#     theta = pi*2**(n-r+1)
#
#     gates = []
#     gates += [Qgate('#', " performs Ry(2^({})*pi/2^({}-1)) = Ry(pi/{}) on {}, controlled by {}"
#                          .format(n, r, 2**(r-n-1), qnb, qna))]
#     gates += [Qgate('ry', qnb, theta/4)]
#     gates += [Qgate('cx', qna, qnb)]
#     gates += [Qgate('ry', qnb, -theta/2)]
#     gates += [Qgate('cx', qna, qnb)]
#     gates += [Qgate('ry', qnb, theta/4)]
#
#     cRysubroutine = Qsubroutine(name="cRy", gates=gates)
#
#     return cRysubroutine


def eigenvalue_inversion(qnl, qnm, qnc, sign=1, remove_global_shift=True):
    """outputs the circuit part designed to invert the eigenvalues"""

    def crzz(qnm, sign=1, remove_global_shift=remove_global_shift):
        gates = []
        grgates = []
        gates += [Qgate('#', ' Performs cRzz({}*pi/2^(m-i)) on L, controlled by the i-th qubit of M, for all i in M'.format(sign*2))]
        if remove_global_shift:
            theta = sign*2*pi/(2**len(qnm))
            grgates += [Qgate('h', qnm[0])]
            grgates += [Qgate('rz', qnm[0], theta)]
            grgates += [Qgate('h', qnm[0])]
        if sign == 1:
            gates += grgates
        for i in range(len(qnm)):
            # applies a cRzz(2*pi/2^(m-i)) to the last L qubit, controlled by the i-th M qubit
            # [= Rz(2*pi/2^(m-i)) on the i-th M qubit, up to global phase]
            gates += [Qgate('rz', qnm[i], sign*2*pi/(2**(len(qnm)-i)))]
        if sign == -1:
            gates += grgates
        gates += [Qgate('display')]
        subroutine = Qsubroutine(name='cRzz', gates=gates)
        return subroutine

    def exp_iH0t2m(qnl, qnm, qnc, t, sign=1):
        """outputs the gates to perform cc-exp(i*H0*t/2^m), which is ccRz(t/2^i) on the i-th qubit of M (for all i)"""
        # Watch out: qnl and qnc should be a specific qubit, while qnm should be all qubits of M!

        gates = []
        gates += [Qgate()]
        gates += [Qgate('#', '## Performs exp(i*H0*t/2^m) with t = {}, on M, controlled by {} and {}'.format(sign*t, qnl, qnc))]
        for i in range(len(qnm)):
            gates += ccRz(qna=qnl, qnb=qnc, qnc=qnm[i], theta=sign*t/(2**(i+1)))
        return gates

    def G(qnl, qnm, qnc, l, kl, t, sign=1):
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
            gates += exp_iH0t2m(qnl=qnl, qnm=qnm, qnc=qnc[kc], t=sign*t/(2**(u + v - m)))
        return gates

    def exp_iH0t0(qnl, qnm, qnc, sign=1):
        """outputs the gates to perform exp(i*H0*t0), which is (qnl[kl], qnm, qnc, l, kl, t0) for all kl in L
        with t0=2*pi"""

        t0 = sign*2*pi
        l = len(qnl)

        gates = []
        gates += [Qgate()]
        gates += [Qgate('#', '###### Performs exp(i*H0*t0) with t0 = {}*pi'.format(sign*2))]
        for kl in range(len(qnc)):
            gates += G(qnl[kl], qnm, qnc, l, kl, t0)
        return gates

    crzzsubroutine = crzz(qnm=qnm, sign=sign, remove_global_shift=remove_global_shift)

    expiH0t0gates = exp_iH0t0(qnl=qnl, qnm=qnm, qnc=qnc, sign=sign)
    expiH0t0subroutine = Qsubroutine(name='exp_iH0t0', gates=expiH0t0gates)

    if sign == 1:
        subroutines = [crzzsubroutine, expiH0t0subroutine]
    elif sign == -1:
        subroutines = [expiH0t0subroutine, crzzsubroutine]
    else:
        raise ValueError('Variable "sign" must be either +1 or -1')

    return subroutines


class EigenvalueInversion_Circuit(Qfunction):

    def __init__(self, n=4, x=2, remove_global_shift=False, test_undo=False, save_value=False):
        # n: size of C
        # x: which qubit of C should be flipped

        qnl = buildnames(n, qubitnames="l")
        qnm = buildnames(n + 4, qubitnames="m")
        qnc = buildnames(n, qubitnames="c")
        if save_value:
            qnd = buildnames(n, qubitnames="d")
            qn = qnl + qnm + qnc + qnd
        else:
            qnd = []
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

        eigenvalueinversionsubroutine = eigenvalue_inversion(
            qnl=qnl,
            qnm=qnm,
            qnc=qnc,
            sign=1,
            remove_global_shift=remove_global_shift)

        measuregates = [Qgate('display'), Qgate('measure')]
        measuresubroutine = Qsubroutine(name='measure', gates=measuregates)

        if test_undo or save_value:
            uneigenvalueinversionsubroutine = eigenvalue_inversion(
                qnl=qnl,
                qnm=qnm,
                qnc=qnc,
                sign=-1,
                remove_global_shift=remove_global_shift)
            uninitgates = []
            for i in range(len(qnl)):
                uninitgates += [Qgate('h', qnl[i])]
            for i in range(len(qnm)):
                uninitgates += [Qgate('h', qnm[i])]
            uninitsubroutine = Qsubroutine(name="uninit", gates=uninitgates)

            if save_value:
                savegates = []
                savegates += [Qgate('#', ' saves the found inverted eigenvalue(s) to the D register')]
                for i in range(len(qnd)):
                    savegates += [Qgate('cx', qnl[i], qnd[i])]
                savesubroutine = Qsubroutine(name='save', gates=savegates)

                subroutines = \
                    [initsubroutine] + \
                    eigenvalueinversionsubroutine + \
                    [savesubroutine] + \
                    uneigenvalueinversionsubroutine + \
                    [uninitsubroutine] + \
                    [measuresubroutine]
            else:
                subroutines = \
                    [initsubroutine] + \
                    eigenvalueinversionsubroutine + \
                    uneigenvalueinversionsubroutine + \
                    [uninitsubroutine] + \
                    [measuresubroutine]
        else:
            subroutines = \
                [initsubroutine] + \
                eigenvalueinversionsubroutine + \
                [measuresubroutine]

        super().__init__(name='Eigenvalue inversion test', qubits=len(qn), subroutines=subroutines)


if __name__ == "__main__":
    circ = EigenvalueInversion_Circuit(n=4, x=2)
    print(circ)

    print('\n\n\n\n\n')

    test_ccrz = Qfunction(
        name='ccRz test',
        qubits=3,
        subroutines=
            Qsubroutine(
                name='test',
                gates=
                    [Qgate(), Qgate('#', ' init'), Qgate('x', 'q0'), Qgate('x', 'q1'), Qgate('display')] +
                    ccRz(
                        qna='q0',
                        qnb='q1',
                        qnc='q2',
                        theta=4.0,
                        different_comment_names=['a', 'b', 'c']) +
                    [Qgate(), Qgate('#', ' meas'), Qgate('measure'), Qgate('display')]))
    print(test_ccrz)
