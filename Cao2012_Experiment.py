from cQASM import *
from QFT import QFT, iQFT, REVERSE
from math import pi


class Cao2012Experiment(Qfunction):

    def __init__(self, r=5):

        def expA(qubitnamea, qubitnameb, qubitnamec, sign=1, n=0):
            """outputs exp(sign*i*A*t0*2^n/16) = exp(sign*2*pi*i*A*2^(n-4))"""

            qna = qubitnamea
            qnb = qubitnameb
            qnc = qubitnamec

            # ccZ(a,b,c)
            ccz_a_b_c_gates = []
            ccz_a_b_c_gates += [Qgate()]
            ccz_a_b_c_gates += [Qgate('#', ' ccZ(a,b,c)')]
            ccz_a_b_c_gates += [Qgate('tdag', qnb)]
            ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
            ccz_a_b_c_gates += [Qgate('t', qnc)]
            ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
            ccz_a_b_c_gates += [Qgate('cx', qna, qnb)]
            ccz_a_b_c_gates += [Qgate('t', qnb)]
            ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
            ccz_a_b_c_gates += [Qgate('tdag', qnc)]
            ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
            ccz_a_b_c_gates += [Qgate('cx', qna, qnb)]
            ccz_a_b_c_gates += [Qgate('tdag', qna)]
            ccz_a_b_c_gates += [Qgate('tdag', qnc)]
            ccz_a_b_c_gates += [Qgate('cx', qna, qnc)]
            ccz_a_b_c_gates += [Qgate('t', qnc)]
            ccz_a_b_c_gates += [Qgate('cx', qna, qnc)]

            # cRx(a,c,sign*0.196)
            crx_a_c_020_gates = []
            crx_a_c_020_gates += [Qgate()]
            crx_a_c_020_gates += [Qgate('#', ' cRx(a,c,{})'.format(-sign*0.196))]
            crx_a_c_020_gates += [Qgate('rx', qnc, -sign*0.196/4)]
            crx_a_c_020_gates += [Qgate('cz', qna, qnc)]
            crx_a_c_020_gates += [Qgate('rx', qnc, sign*0.196/2)]
            crx_a_c_020_gates += [Qgate('cz', qna, qnc)]
            crx_a_c_020_gates += [Qgate('rx', qnc, -sign*0.196/4)]

            # cVdag(a,c)
            cvdag_a_c_gates = []
            cvdag_a_c_gates += [Qgate()]
            cvdag_a_c_gates += [Qgate('#', ' cVdag(a,c)')]
            # cH(a,c)
            cvdag_a_c_gates += [Qgate('t', qna)]
            cvdag_a_c_gates += [Qgate('x', qna)]
            cvdag_a_c_gates += [Qgate('tdag', qna)]
            cvdag_a_c_gates += [Qgate('x', qna)]
            cvdag_a_c_gates += [Qgate('h', qnc)]
            cvdag_a_c_gates += [Qgate('tdag', qnc)]
            cvdag_a_c_gates += [Qgate('tdag', qnc)]
            cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
            cvdag_a_c_gates += [Qgate('h', qnc)]
            cvdag_a_c_gates += [Qgate('t', qnc)]
            cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
            cvdag_a_c_gates += [Qgate('t', qnc)]
            cvdag_a_c_gates += [Qgate('h', qnc)]
            cvdag_a_c_gates += [Qgate('s', qnc)]
            cvdag_a_c_gates += [Qgate('x', qnc)]
            # cZ(a,c)
            cvdag_a_c_gates += [Qgate('cz', qna, qnc)]
            # cS(a,c)
            cvdag_a_c_gates += [Qgate('t', qna)]
            cvdag_a_c_gates += [Qgate('t', qnc)]
            cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
            cvdag_a_c_gates += [Qgate('tdag', qnc)]
            cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
            # cH(a,c)
            cvdag_a_c_gates += [Qgate('t', qna)]
            cvdag_a_c_gates += [Qgate('x', qna)]
            cvdag_a_c_gates += [Qgate('tdag', qna)]
            cvdag_a_c_gates += [Qgate('x', qna)]
            cvdag_a_c_gates += [Qgate('h', qnc)]
            cvdag_a_c_gates += [Qgate('tdag', qnc)]
            cvdag_a_c_gates += [Qgate('tdag', qnc)]
            cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
            cvdag_a_c_gates += [Qgate('h', qnc)]
            cvdag_a_c_gates += [Qgate('t', qnc)]
            cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
            cvdag_a_c_gates += [Qgate('t', qnc)]
            cvdag_a_c_gates += [Qgate('h', qnc)]
            cvdag_a_c_gates += [Qgate('s', qnc)]
            cvdag_a_c_gates += [Qgate('x', qnc)]

            # cV(a,c)
            cv_a_c_gates = []
            cv_a_c_gates += [Qgate()]
            cv_a_c_gates += [Qgate('#', ' cV(a,c)')]
            # cH(a,c)
            cv_a_c_gates += [Qgate('t', qna)]
            cv_a_c_gates += [Qgate('x', qna)]
            cv_a_c_gates += [Qgate('tdag', qna)]
            cv_a_c_gates += [Qgate('x', qna)]
            cv_a_c_gates += [Qgate('h', qnc)]
            cv_a_c_gates += [Qgate('tdag', qnc)]
            cv_a_c_gates += [Qgate('tdag', qnc)]
            cv_a_c_gates += [Qgate('cx', qna, qnc)]
            cv_a_c_gates += [Qgate('h', qnc)]
            cv_a_c_gates += [Qgate('t', qnc)]
            cv_a_c_gates += [Qgate('cx', qna, qnc)]
            cv_a_c_gates += [Qgate('t', qnc)]
            cv_a_c_gates += [Qgate('h', qnc)]
            cv_a_c_gates += [Qgate('s', qnc)]
            cv_a_c_gates += [Qgate('x', qnc)]
            # cS(a,c)
            cv_a_c_gates += [Qgate('t', qna)]
            cv_a_c_gates += [Qgate('t', qnc)]
            cv_a_c_gates += [Qgate('cx', qna, qnc)]
            cv_a_c_gates += [Qgate('tdag', qnc)]
            cv_a_c_gates += [Qgate('cx', qna, qnc)]
            # cH(a,c)
            cv_a_c_gates += [Qgate('t', qna)]
            cv_a_c_gates += [Qgate('x', qna)]
            cv_a_c_gates += [Qgate('tdag', qna)]
            cv_a_c_gates += [Qgate('x', qna)]
            cv_a_c_gates += [Qgate('h', qnc)]
            cv_a_c_gates += [Qgate('tdag', qnc)]
            cv_a_c_gates += [Qgate('tdag', qnc)]
            cv_a_c_gates += [Qgate('cx', qna, qnc)]
            cv_a_c_gates += [Qgate('h', qnc)]
            cv_a_c_gates += [Qgate('t', qnc)]
            cv_a_c_gates += [Qgate('cx', qna, qnc)]
            cv_a_c_gates += [Qgate('t', qnc)]
            cv_a_c_gates += [Qgate('h', qnc)]
            cv_a_c_gates += [Qgate('s', qnc)]
            cv_a_c_gates += [Qgate('x', qnc)]

            # Rz(a,sign*0.375)
            rz_a_038_gates = []
            rz_a_038_gates += [Qgate()]
            rz_a_038_gates += [Qgate('#', ' Rz(a,{})'.format(sign*0.375))]
            rz_a_038_gates += [Qgate('rz', qna, sign*0.375)]

            # cRx(a,b,-sign*0.982)
            crx_a_b_098_gates = []
            crx_a_b_098_gates += [Qgate()]
            crx_a_b_098_gates += [Qgate('#', ' cRx(a,b,{})'.format(-sign*0.982))]
            crx_a_b_098_gates += [Qgate('rx', qnb, -sign*0.982/4)]
            crx_a_b_098_gates += [Qgate('cz', qna, qnb)]
            crx_a_b_098_gates += [Qgate('rx', qnb, sign*0.982/2)]
            crx_a_b_098_gates += [Qgate('cz', qna, qnb)]
            crx_a_b_098_gates += [Qgate('rx', qnb, -sign*0.982/4)]

            # Rz(a,sign*1.883)
            rz_a_188_gates = []
            rz_a_188_gates += [Qgate()]
            rz_a_188_gates += [Qgate('#', ' Rz(a,{})'.format(sign*1.883))]
            rz_a_188_gates += [Qgate('rz', qna, sign*1.883)]

            # Toffoli(a,b,c)
            toffoli_a_b_c_gates = []
            toffoli_a_b_c_gates += [Qgate()]
            toffoli_a_b_c_gates += [Qgate('#', ' Toffoli(a,b,c)')]
            toffoli_a_b_c_gates += [Qgate('toffoli', qna, qnb, qnc)]

            # cRx(a,b,-sign*0.589)
            crx_a_b_059_gates = []
            crx_a_b_059_gates += [Qgate()]
            crx_a_b_059_gates += [Qgate('#', ' cRx(a,b,{})'.format(-sign*0.589))]
            crx_a_b_059_gates += [Qgate('rx', qnb, -sign*0.589/4)]
            crx_a_b_059_gates += [Qgate('cz', qna, qnb)]
            crx_a_b_059_gates += [Qgate('rx', qnb, sign*0.589/2)]
            crx_a_b_059_gates += [Qgate('cz', qna, qnb)]
            crx_a_b_059_gates += [Qgate('rx', qnb, -sign*0.589/4)]

            # # Toffoli(a,b,c)
            # toffoli_a_b_c_gates = []
            # toffoli_a_b_c_gates += [Qgate()]
            # toffoli_a_b_c_gates += [Qgate('#', ' Toffoli(a,b,c)')]
            # toffoli_a_b_c_gates += [Qgate('toffoli', qna, qnb, qnc)]

            # cZ(a,c)
            cz_a_c_gates = []
            cz_a_c_gates += [Qgate()]
            cz_a_c_gates += [Qgate('#', ' cZ(a,c)')]
            cz_a_c_gates += [Qgate('cz', qna, qnc)]

            # cX(a,b)
            cx_a_b_gates = []
            cx_a_b_gates += [Qgate()]
            cx_a_b_gates += [Qgate('#', ' cX(a,b)')]
            cx_a_b_gates += [Qgate('cx', qna, qnb)]

            gates = []

            if sign == 1:
                gates += ccz_a_b_c_gates + \
                         crx_a_c_020_gates + \
                         cvdag_a_c_gates + \
                         rz_a_038_gates + \
                         crx_a_b_098_gates + \
                         rz_a_188_gates + \
                         toffoli_a_b_c_gates + \
                         crx_a_b_059_gates + \
                         toffoli_a_b_c_gates + \
                         cz_a_c_gates + \
                         cx_a_b_gates + \
                         ccz_a_b_c_gates + \
                         cx_a_b_gates
                
            elif sign == -1:
                gates += cx_a_b_gates + \
                         ccz_a_b_c_gates + \
                         cx_a_b_gates + \
                         cz_a_c_gates + \
                         toffoli_a_b_c_gates + \
                         crx_a_b_059_gates + \
                         toffoli_a_b_c_gates + \
                         rz_a_188_gates + \
                         crx_a_b_098_gates + \
                         rz_a_038_gates + \
                         cv_a_c_gates + \
                         crx_a_c_020_gates + \
                         ccz_a_b_c_gates
                
            else:
                raise ValueError("variable 'sign' must be either +1 or -1")

            for i in range(n):
                gates += [Qgate(), Qgate('#', '####')] + gates

            gates.insert(0,
                         Qgate('#', " performs exp(sign*i*A*t0*2^(n)/16) = exp({}*pi*i*A/{}) "
                                    "on {} and {}, controlled by {}"
                                    .format(2*sign, 2**(4-n), qnb, qnc, qna)))

            expAsubroutine = Qsubroutine(name="expA", gates=gates)

            return expAsubroutine

        def cRy(qubitnamea, qubitnameb, n, r):
            """"outputs controlled-Ry(2^(n)*pi/2^(r-1)) on qubit b, with qubit a as control"""

            qna = qubitnamea
            qnb = qubitnameb

            theta = pi*2**(n-r+1)

            gates = []
            gates += [Qgate('#', " performs Ry(-2^({})*pi/2^({}-1)) = Ry(-pi/{}) on {}, controlled by {}"
                                 .format(n, r, 2**(r-n-1), qnb, qna))]
            gates += [Qgate('ry', qnb, -theta/4)]
            gates += [Qgate('cx', qna, qnb)]
            gates += [Qgate('ry', qnb, theta/2)]
            gates += [Qgate('cx', qna, qnb)]
            gates += [Qgate('ry', qnb, -theta/4)]

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

        # reversesubroutine = REVERSE(n=4, qubitnames=qn[1:5])
        iqftsubroutine = iQFT(n=4, qubitnames=qn[1:5])

        swapsubroutine = Qsubroutine(name="swap", gates=[Qgate("swap", qn[2], qn[4])])

        cRysubroutines = []
        for i in range(4):
            cRysubroutines += [cRy(qubitnamea=qn[i+1], qubitnameb=qn[0], n=3-i, r=r)]

        qftsubroutine = QFT(n=4, qubitnames=qn[1:5])

        unexpAsubroutines = []
        for i in reversed(range(4)):
            unexpAsubroutines += [expA(qubitnamea=qn[i+1], qubitnameb=qn[-2], qubitnamec=qn[-1], sign=1, n=i)]

        uninitgates = []
        for i in range(1, 5):
            uninitgates += [Qgate("h", qn[i])]
        uninitsubroutine = Qsubroutine(name="uninit", gates=uninitgates)

        #resultgates = [Qgate("measure"), Qgate("display"), Qgate("measure"), Qgate("display")]
        #resultgates[0].options = ["q0"]
        resultgates = [Qgate("display"), Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = []
        subroutines += [initsubroutine]
        subroutines += expAsubroutines
        subroutines += [iqftsubroutine]
        # subroutines += [reversesubroutine]
        subroutines += [swapsubroutine]
        subroutines += cRysubroutines
        subroutines += [swapsubroutine]
        # subroutines += [reversesubroutine]
        subroutines += [qftsubroutine]
        subroutines += unexpAsubroutines
        subroutines += [uninitsubroutine]
        subroutines += [resultsubroutine]

        super().__init__(name="Cao2012 Experiment", qubits=qubits, subroutines=subroutines)
