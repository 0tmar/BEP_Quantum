from cQASM import *
from QFT import QFT, iQFT
from math import pi


class Cao2012Experiment(Qfunction):

    def __init__(self, r=5):

        def expA(qubitnamea, qubitnameb, qubitnamec, sign=1, n=0):
            """outputs exp(sign*i*A*t0*2^n/16) = exp(sign*2*pi*i*A*2^(n-4))"""

            qna = qubitnamea
            qnb = qubitnameb
            qnc = qubitnamec

            gates = []

            # ccZ(a,b,c)
            gates += [Qgate()]
            gates += [Qgate('#', ' ccZ(a,b,c)')]
            gates += [Qgate('tdag', qnb)]
            gates += [Qgate('cx', qnb, qnc)]
            gates += [Qgate('t', qnc)]
            gates += [Qgate('cx', qnb, qnc)]
            gates += [Qgate('cx', qna, qnb)]
            gates += [Qgate('t', qnb)]
            gates += [Qgate('cx', qnb, qnc)]
            gates += [Qgate('tdag', qnc)]
            gates += [Qgate('cx', qnb, qnc)]
            gates += [Qgate('cx', qna, qnb)]
            gates += [Qgate('tdag', qna)]
            gates += [Qgate('tdag', qnc)]
            gates += [Qgate('cx', qna, qnc)]
            gates += [Qgate('t', qnc)]
            gates += [Qgate('cx', qna, qnc)]

            # cRx(a,c,0.2)
            gates += [Qgate()]
            gates += [Qgate('#', ' cRx(a,c,0.2)')]
            gates += [Qgate('rx', qnc, 0.2/4)]
            gates += [Qgate('cz', qna, qnc)]
            gates += [Qgate('rx', qnc, -0.2/2)]
            gates += [Qgate('cz', qna, qnc)]
            gates += [Qgate('rx', qnc, 0.2/4)]

            # cVdag(a,c)
            gates += [Qgate()]
            gates += [Qgate('#', ' cVdag(a,c)')]
            #   cH(a,c)
            gates += [Qgate('s', qna)]
            gates += [Qgate('h', qnc)]
            gates += [Qgate('tdag', qnc)]
            gates += [Qgate('tdag', qnc)]
            gates += [Qgate('cx', qna, qnc)]
            gates += [Qgate('h', qnc)]
            gates += [Qgate('t', qnc)]
            gates += [Qgate('cx', qna, qnc)]
            gates += [Qgate('t', qnc)]
            gates += [Qgate('h', qnc)]
            gates += [Qgate('s', qnc)]
            gates += [Qgate('x', qnc)]
            #   cZ(a,c)
            gates += [Qgate('cz', qna, qnc)]
            #   cS(a,c)
            gates += [Qgate('tdag', qna)]
            gates += [Qgate('tdag', qnc)]
            gates += [Qgate('cx', qna, qnc)]
            gates += [Qgate('t', qnc)]
            gates += [Qgate('cx', qna, qnc)]
            #   cH(a,c)
            gates += [Qgate('s', qna)]
            gates += [Qgate('h', qnc)]
            gates += [Qgate('tdag', qnc)]
            gates += [Qgate('tdag', qnc)]
            gates += [Qgate('cx', qna, qnc)]
            gates += [Qgate('h', qnc)]
            gates += [Qgate('t', qnc)]
            gates += [Qgate('cx', qna, qnc)]
            gates += [Qgate('t', qnc)]
            gates += [Qgate('h', qnc)]
            gates += [Qgate('s', qnc)]
            gates += [Qgate('x', qnc)]

            # Rz(a,0.38)
            gates += [Qgate()]
            gates += [Qgate('#', ' Rz(a,0.38)')]
            gates += [Qgate('rz', qna, 0.38)]

            # cRx(a,b,0.98)
            gates += [Qgate()]
            gates += [Qgate('#', ' cRx(a,b,0.98)')]
            gates += [Qgate('rx', qnc, 0.98/4)]
            gates += [Qgate('cz', qna, qnc)]
            gates += [Qgate('rx', qnc, -0.98/2)]
            gates += [Qgate('cz', qna, qnc)]
            gates += [Qgate('rx', qnc, 0.98/4)]

            # Rz(a,1.88)
            gates += [Qgate()]
            gates += [Qgate('#', ' Rz(a,1.88)')]
            gates += [Qgate('rz', qna, 1.88)]

            # Toffoli(a,b,c)
            gates += [Qgate()]
            gates += [Qgate('#', ' Toffoli(a,b,c)')]
            gates += [Qgate('toffoli', qna, qnb, qnc)]

            # cRx(a,b,0.59)
            gates += [Qgate()]
            gates += [Qgate('#', ' cRx(a,b,0.59)')]
            gates += [Qgate('rx', qnc, 0.59/4)]
            gates += [Qgate('cz', qna, qnc)]
            gates += [Qgate('rx', qnc, -0.59/2)]
            gates += [Qgate('cz', qna, qnc)]
            gates += [Qgate('rx', qnc, 0.59/4)]

            # Toffoli(a,b,c)
            gates += [Qgate()]
            gates += [Qgate('#', ' Toffoli(a,b,c)')]
            gates += [Qgate('toffoli', qna, qnb, qnc)]

            # cZ(a,c)
            gates += [Qgate()]
            gates += [Qgate('#', ' cZ(a,c)')]
            gates += [Qgate('cz', qna, qnc)]

            for i in range(n):
                gates += [Qgate(), Qgate('#', '####')] + gates

            gates.insert(0,
                         Qgate('#', " performs exp(sign*i*A*t0*2^(n)/16) = exp({}*pi*i*A/{}) "
                                    "on {} and {}, controlled by {}"
                                    .format(2*sign, 2**(4-n), qnb, qnc, qna)))

            expAsr = Qsubroutine(name="expA", gates=gates)

            return expAsr

        def cRy(qubitnamea, qubitnameb, n, r):
            """"outputs controlled-Ry(2^(n)*pi/2^(r-1)) on qubit b, with qubit a as control"""

            qna = qubitnamea
            qnb = qubitnameb

            theta = pi*2**(n-r+1)

            gates = []
            gates += [Qgate('#', " performs Ry(2^({})*pi/2^({}-1)) = Ry(pi/2^{}) on {}, controlled by {}"
                                 .format(n, r, 2**(r-n-1), qnb, qna))]
            gates += [Qgate('ry', qnb, theta/4)]
            gates += [Qgate('cx', qna, qnb)]
            gates += [Qgate('ry', qnb, -theta/2)]
            gates += [Qgate('cx', qna, qnb)]
            gates += [Qgate('ry', qnb, theta/4)]

            cRysr = Qsubroutine(name="cRy", gates=gates)

            return cRysr

        qubits = 7

        qn = buildnames(n=7, qubitnames="q")

        initgates = []
        for i in range(1, 7):
            initgates += [Qgate("h", qn[i])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        expAsubroutines = []
        for i in range(4):
            expAsubroutines += [expA(qubitnamea=qn[i+1], qubitnameb=qn[-2], qubitnamec=qn[-1], sign=-1, n=i)]

        iqftsubroutine = iQFT(n=4, qubitnames=qn[1:5])

        swapsubroutine = Qsubroutine(name="swap", gates=[Qgate("swap", qn[2], qn[4])])

        cRysubroutines = []
        for i in range(4):
            cRysubroutines += [cRy(qubitnamea=qn[i+1], qubitnameb=qn[0], n=3-i, r=r)]

        qftsubroutine = QFT(n=4, qubitnames=qn[1:5])

        unexpAsubroutines = []
        for i in reversed(range(4)):
            unexpAsubroutines += [expA(qubitnamea=qn[i + 1], qubitnameb=qn[-2], qubitnamec=qn[-1], sign=1, n=i)]

        uninitgates = []
        for i in range(1, 5):
            uninitgates += [Qgate("h", qn[i])]
        uninitsubroutine = Qsubroutine(name="uninit", gates=uninitgates)

        resultgates = [Qgate("display"), Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = []
        subroutines += [initsubroutine]
        subroutines += expAsubroutines
        subroutines += [iqftsubroutine]
        subroutines += [swapsubroutine]
        subroutines += cRysubroutines
        subroutines += [swapsubroutine]
        subroutines += [qftsubroutine]
        subroutines += unexpAsubroutines
        subroutines += [uninitsubroutine]
        subroutines += [resultsubroutine]

        super().__init__(name="Cao2012 Experiment", qubits=qubits, subroutines=subroutines)
