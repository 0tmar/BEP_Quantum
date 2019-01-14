from cQASM import *
from cQASM_Compound_Gates import *
from QFT import QFT, iQFT, REVERSE
from math import pi


def expA(qubitnamea, qubitnameb, qubitnamec, qubitnamed=None, sign=1, n=0, noglobalrotation=False):
    """outputs exp(sign*i*A*t0*2^n/16) = exp(sign*2*pi*i*A*2^(n-4))"""

    qna = qubitnamea
    qnb = qubitnameb
    qnc = qubitnamec

    # ccZ(a,b,c)
    ccz_a_b_c_gates = ccZ(qna=qna, qnb=qnb, qnc=qnc, different_comment_names=['a', 'b', 'c'])

    # cRx(a,c,-sign*(2**n)*0.196)
    crx_a_c_020_gates = cRx(qna=qna, qnb=qnc, theta=-sign*(2**n)*0.196, different_comment_names=['a', 'c'])

    # cVdag(a,c)
    cvdag_a_c_gates = cVdag(qna=qna, qnb=qnc, different_comment_names=['a', 'c'])

    # cV(a,c)
    cv_a_c_gates = cV(qna=qna, qnb=qnc, different_comment_names=['a', 'c'])

    # cX(a,c)
    cx_a_c_gates = []
    cx_a_c_gates += [Qgate()]
    cx_a_c_gates += [Qgate('#', ' cX(a,c)')]
    cx_a_c_gates += [Qgate('cx', qna, qnc)]

    # Rz(a,sign*(2**n)*0.375)
    rz_a_038_gates = []
    rz_a_038_gates += [Qgate()]
    rz_a_038_gates += [Qgate('#', ' Rz(a,{})'.format(sign * (2 ** n) * 0.375))]
    rz_a_038_gates += [Qgate('rz', qna, sign * (2 ** n) * 0.375)]

    # cRx(a,b,-sign*(2**n)*0.982)
    crx_a_b_098_gates = cRx(qna=qna, qnb=qnb, theta=-sign*(2**n)*0.982, different_comment_names=['a', 'b'])

    # Rz(a,sign*(2**n)*1.883)
    rz_a_188_gates = []
    rz_a_188_gates += [Qgate()]
    rz_a_188_gates += [Qgate('#', ' Rz(a,{})'.format(sign * (2 ** n) * 1.883))]
    rz_a_188_gates += [Qgate('rz', qna, sign * (2 ** n) * 1.883)]

    # Toffoli(a,b,c)
    toffoli_a_b_c_gates = []
    toffoli_a_b_c_gates += [Qgate()]
    toffoli_a_b_c_gates += [Qgate('#', ' Toffoli(a,b,c)')]
    toffoli_a_b_c_gates += [Qgate('toffoli', qna, qnb, qnc)]

    # cRx(a,b,-sign*(2**n)*0.589)
    crx_a_b_059_gates = cRx(qna=qna, qnb=qnb, theta=-sign*(2**n)*0.589, different_comment_names=['a', 'b'])

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

    # # ccZ(a,b,c)
    # ccz_a_b_c_gates = []
    # ccz_a_b_c_gates += [Qgate()]
    # ccz_a_b_c_gates += [Qgate('#', ' ccZ(a,b,c)')]
    # ccz_a_b_c_gates += [Qgate('tdag', qnb)]
    # ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
    # ccz_a_b_c_gates += [Qgate('t', qnc)]
    # ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
    # ccz_a_b_c_gates += [Qgate('cx', qna, qnb)]
    # ccz_a_b_c_gates += [Qgate('t', qnb)]
    # ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
    # ccz_a_b_c_gates += [Qgate('tdag', qnc)]
    # ccz_a_b_c_gates += [Qgate('cx', qnb, qnc)]
    # ccz_a_b_c_gates += [Qgate('cx', qna, qnb)]
    # ccz_a_b_c_gates += [Qgate('tdag', qna)]
    # ccz_a_b_c_gates += [Qgate('tdag', qnc)]
    # ccz_a_b_c_gates += [Qgate('cx', qna, qnc)]
    # ccz_a_b_c_gates += [Qgate('t', qnc)]
    # ccz_a_b_c_gates += [Qgate('cx', qna, qnc)]
    #
    # # cRx(a,c,-sign*(2**n)*0.196)
    # crx_a_c_020_gates = []
    # crx_a_c_020_gates += [Qgate()]
    # crx_a_c_020_gates += [Qgate('#', ' cRx(a,c,{})'.format(-sign * (2 ** n) * 0.196))]
    # crx_a_c_020_gates += [Qgate('rx', qnc, -sign * (2 ** n) * 0.196 / 4)]
    # crx_a_c_020_gates += [Qgate('cz', qna, qnc)]
    # crx_a_c_020_gates += [Qgate('rx', qnc, sign * (2 ** n) * 0.196 / 2)]
    # crx_a_c_020_gates += [Qgate('cz', qna, qnc)]
    # crx_a_c_020_gates += [Qgate('rx', qnc, -sign * (2 ** n) * 0.196 / 4)]
    #
    # # cVdag(a,c)
    # cvdag_a_c_gates = []
    # cvdag_a_c_gates += [Qgate()]
    # cvdag_a_c_gates += [Qgate('#', ' cVdag(a,c)')]
    # # cH(a,c)
    # cvdag_a_c_gates += [Qgate('t', qna)]
    # cvdag_a_c_gates += [Qgate('x', qna)]
    # cvdag_a_c_gates += [Qgate('tdag', qna)]
    # cvdag_a_c_gates += [Qgate('x', qna)]
    # cvdag_a_c_gates += [Qgate('h', qnc)]
    # cvdag_a_c_gates += [Qgate('tdag', qnc)]
    # cvdag_a_c_gates += [Qgate('tdag', qnc)]
    # cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
    # cvdag_a_c_gates += [Qgate('h', qnc)]
    # cvdag_a_c_gates += [Qgate('t', qnc)]
    # cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
    # cvdag_a_c_gates += [Qgate('t', qnc)]
    # cvdag_a_c_gates += [Qgate('h', qnc)]
    # cvdag_a_c_gates += [Qgate('s', qnc)]
    # cvdag_a_c_gates += [Qgate('x', qnc)]
    # # cSdag(a,c)
    # cvdag_a_c_gates += [Qgate('tdag', qna)]
    # cvdag_a_c_gates += [Qgate('tdag', qnc)]
    # cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
    # cvdag_a_c_gates += [Qgate('t', qnc)]
    # cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
    # # cH(a,c)
    # cvdag_a_c_gates += [Qgate('t', qna)]
    # cvdag_a_c_gates += [Qgate('x', qna)]
    # cvdag_a_c_gates += [Qgate('tdag', qna)]
    # cvdag_a_c_gates += [Qgate('x', qna)]
    # cvdag_a_c_gates += [Qgate('h', qnc)]
    # cvdag_a_c_gates += [Qgate('tdag', qnc)]
    # cvdag_a_c_gates += [Qgate('tdag', qnc)]
    # cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
    # cvdag_a_c_gates += [Qgate('h', qnc)]
    # cvdag_a_c_gates += [Qgate('t', qnc)]
    # cvdag_a_c_gates += [Qgate('cx', qna, qnc)]
    # cvdag_a_c_gates += [Qgate('t', qnc)]
    # cvdag_a_c_gates += [Qgate('h', qnc)]
    # cvdag_a_c_gates += [Qgate('s', qnc)]
    # cvdag_a_c_gates += [Qgate('x', qnc)]
    #
    # # cV(a,c)
    # cv_a_c_gates = []
    # cv_a_c_gates += [Qgate()]
    # cv_a_c_gates += [Qgate('#', ' cV(a,c)')]
    # # cH(a,c)
    # cv_a_c_gates += [Qgate('t', qna)]
    # cv_a_c_gates += [Qgate('x', qna)]
    # cv_a_c_gates += [Qgate('tdag', qna)]
    # cv_a_c_gates += [Qgate('x', qna)]
    # cv_a_c_gates += [Qgate('h', qnc)]
    # cv_a_c_gates += [Qgate('tdag', qnc)]
    # cv_a_c_gates += [Qgate('tdag', qnc)]
    # cv_a_c_gates += [Qgate('cx', qna, qnc)]
    # cv_a_c_gates += [Qgate('h', qnc)]
    # cv_a_c_gates += [Qgate('t', qnc)]
    # cv_a_c_gates += [Qgate('cx', qna, qnc)]
    # cv_a_c_gates += [Qgate('t', qnc)]
    # cv_a_c_gates += [Qgate('h', qnc)]
    # cv_a_c_gates += [Qgate('s', qnc)]
    # cv_a_c_gates += [Qgate('x', qnc)]
    # # cS(a,c)
    # cv_a_c_gates += [Qgate('t', qna)]
    # cv_a_c_gates += [Qgate('t', qnc)]
    # cv_a_c_gates += [Qgate('cx', qna, qnc)]
    # cv_a_c_gates += [Qgate('tdag', qnc)]
    # cv_a_c_gates += [Qgate('cx', qna, qnc)]
    # # cH(a,c)
    # cv_a_c_gates += [Qgate('t', qna)]
    # cv_a_c_gates += [Qgate('x', qna)]
    # cv_a_c_gates += [Qgate('tdag', qna)]
    # cv_a_c_gates += [Qgate('x', qna)]
    # cv_a_c_gates += [Qgate('h', qnc)]
    # cv_a_c_gates += [Qgate('tdag', qnc)]
    # cv_a_c_gates += [Qgate('tdag', qnc)]
    # cv_a_c_gates += [Qgate('cx', qna, qnc)]
    # cv_a_c_gates += [Qgate('h', qnc)]
    # cv_a_c_gates += [Qgate('t', qnc)]
    # cv_a_c_gates += [Qgate('cx', qna, qnc)]
    # cv_a_c_gates += [Qgate('t', qnc)]
    # cv_a_c_gates += [Qgate('h', qnc)]
    # cv_a_c_gates += [Qgate('s', qnc)]
    # cv_a_c_gates += [Qgate('x', qnc)]
    #
    # # cX(a,c)
    # cx_a_c_gates = []
    # cx_a_c_gates += [Qgate()]
    # cx_a_c_gates += [Qgate('#', ' cX(a,c)')]
    # cx_a_c_gates += [Qgate('cx', qna, qnc)]
    #
    # # Rz(a,sign*(2**n)*0.375)
    # rz_a_038_gates = []
    # rz_a_038_gates += [Qgate()]
    # rz_a_038_gates += [Qgate('#', ' Rz(a,{})'.format(sign * (2 ** n) * 0.375))]
    # rz_a_038_gates += [Qgate('rz', qna, sign * (2 ** n) * 0.375)]
    #
    # # cRx(a,b,-sign*(2**n)*0.982)
    # crx_a_b_098_gates = []
    # crx_a_b_098_gates += [Qgate()]
    # crx_a_b_098_gates += [Qgate('#', ' cRx(a,b,{})'.format(-sign * (2 ** n) * 0.982))]
    # crx_a_b_098_gates += [Qgate('rx', qnb, -sign * (2 ** n) * 0.982 / 4)]
    # crx_a_b_098_gates += [Qgate('cz', qna, qnb)]
    # crx_a_b_098_gates += [Qgate('rx', qnb, sign * (2 ** n) * 0.982 / 2)]
    # crx_a_b_098_gates += [Qgate('cz', qna, qnb)]
    # crx_a_b_098_gates += [Qgate('rx', qnb, -sign * (2 ** n) * 0.982 / 4)]
    #
    # # Rz(a,sign*(2**n)*1.883)
    # rz_a_188_gates = []
    # rz_a_188_gates += [Qgate()]
    # rz_a_188_gates += [Qgate('#', ' Rz(a,{})'.format(sign * (2 ** n) * 1.883))]
    # rz_a_188_gates += [Qgate('rz', qna, sign * (2 ** n) * 1.883)]
    #
    # # Toffoli(a,b,c)
    # toffoli_a_b_c_gates = []
    # toffoli_a_b_c_gates += [Qgate()]
    # toffoli_a_b_c_gates += [Qgate('#', ' Toffoli(a,b,c)')]
    # toffoli_a_b_c_gates += [Qgate('toffoli', qna, qnb, qnc)]
    #
    # # cRx(a,b,-sign*(2**n)*0.589)
    # crx_a_b_059_gates = []
    # crx_a_b_059_gates += [Qgate()]
    # crx_a_b_059_gates += [Qgate('#', ' cRx(a,b,{})'.format(-sign * (2 ** n) * 0.589))]
    # crx_a_b_059_gates += [Qgate('rx', qnb, -sign * (2 ** n) * 0.589 / 4)]
    # crx_a_b_059_gates += [Qgate('cz', qna, qnb)]
    # crx_a_b_059_gates += [Qgate('rx', qnb, sign * (2 ** n) * 0.589 / 2)]
    # crx_a_b_059_gates += [Qgate('cz', qna, qnb)]
    # crx_a_b_059_gates += [Qgate('rx', qnb, -sign * (2 ** n) * 0.589 / 4)]
    #
    # # # Toffoli(a,b,c)
    # # toffoli_a_b_c_gates = []
    # # toffoli_a_b_c_gates += [Qgate()]
    # # toffoli_a_b_c_gates += [Qgate('#', ' Toffoli(a,b,c)')]
    # # toffoli_a_b_c_gates += [Qgate('toffoli', qna, qnb, qnc)]
    #
    # # cZ(a,c)
    # cz_a_c_gates = []
    # cz_a_c_gates += [Qgate()]
    # cz_a_c_gates += [Qgate('#', ' cZ(a,c)')]
    # cz_a_c_gates += [Qgate('cz', qna, qnc)]
    #
    # # cX(a,b)
    # cx_a_b_gates = []
    # cx_a_b_gates += [Qgate()]
    # cx_a_b_gates += [Qgate('#', ' cX(a,b)')]
    # cx_a_b_gates += [Qgate('cx', qna, qnb)]

    if noglobalrotation:
        qnd = qubitnamed

        # Rz(d,2*sign*(2**n)*1.129)
        rz_d_113_gates = []
        rz_d_113_gates += [Qgate()]
        rz_d_113_gates += [Qgate('#', ' Rz(d,{})'.format(sign*(2**n)*2*1.129))]
        rz_d_113_gates += [Qgate('rz', qnd, sign*(2**n)*2.258)]
    else:
        rz_d_113_gates = []

    gates = []

    if sign == 1:
        gates += ccz_a_b_c_gates + \
                 crx_a_c_020_gates

        if n == 0:
            gates += cvdag_a_c_gates
        elif n == 1:
            gates += cx_a_c_gates
        elif n == 2:
            pass
        elif n == 3:
            pass

        gates += rz_a_038_gates + \
                 crx_a_b_098_gates + \
                 rz_a_188_gates + \
                 toffoli_a_b_c_gates + \
                 crx_a_b_059_gates + \
                 toffoli_a_b_c_gates + \
                 cz_a_c_gates + \
                 cx_a_b_gates + \
                 ccz_a_b_c_gates + \
                 cx_a_b_gates

        if noglobalrotation:
            gates += rz_d_113_gates

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
                 rz_a_038_gates

        if n == 0:
            gates += cv_a_c_gates
        elif n == 1:
            gates += cx_a_c_gates
        elif n == 2:
            pass
        elif n == 3:
            pass

        gates += crx_a_c_020_gates + \
                 ccz_a_b_c_gates

        if noglobalrotation:
            gates = rz_d_113_gates + gates

    else:
        raise ValueError("variable 'sign' must be either +1 or -1")

    # for i in range(n):
    #    gates += [Qgate(), Qgate('#', '####')] + gates

    gates.insert(0,
                 Qgate('#', " performs exp(sign*i*A*t0*2^(n)/16) = exp({}*pi*i*A/{}) "
                            "on {} and {}, controlled by {}"
                       .format(2 * sign, 2 ** (4 - n), qnb, qnc, qna)))

    expAsubroutine = Qsubroutine(name="expA", gates=gates)

    return expAsubroutine


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

    cRysr = Qsubroutine(name="cRy", gates=gates)

    return cRysr


class Cao2012Experiment(Qfunction):

    def __init__(self, r=5, m=None, n=None):

        qubits = 7

        qn = buildnames(n=7, qubitnames="q")

        initgates = []
        for i in range(1, 5):
            initgates += [Qgate("h", qn[i])]
        if n is None:
            for i in range(5, 7):
                initgates += [Qgate("h", qn[i])]
        if m is not None:
            if m == 0 or m == 1:
                initgates += [Qgate("x", qn[5])]
            else:
                initgates += [Qgate("#", "x {}".format(qn[5]))]
            if m == 0 or m == 2:
                initgates += [Qgate("x", qn[6])]
            else:
                initgates += [Qgate("#", "x {}".format(qn[6]))]
            initgates += [Qgate("cz", qn[5], qn[6])]
            if m == 0 or m == 2:
                initgates += [Qgate("x", qn[6])]
            else:
                initgates += [Qgate("#", "x {}".format(qn[6]))]
            if m == 0 or m == 1:
                initgates += [Qgate("x", qn[5])]
            else:
                initgates += [Qgate("#", "x {}".format(qn[5]))]
        if n is not None:
            if n == 2 or n == 3:
                initgates += [Qgate("x", qn[5])]
            if n == 1 or n == 3:
                initgates += [Qgate("x", qn[6])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        expAsubroutines = []
        for i in range(4):
            expAsubroutines += [expA(qubitnamea=qn[i+1], qubitnameb=qn[5], qubitnamec=qn[6], sign=1, n=i)]

        reversesubroutine = REVERSE(n=4, qubitnames=qn[1:5])

        iqftsubroutine = iQFT(n=4, qubitnames=qn[1:5])
        iqftsubroutine.gates.append(Qgate())
        iqftsubroutine.gates.append(Qgate("display"))
        iqftsubroutine.gates.insert(0, Qgate())
        iqftsubroutine.gates.insert(0, Qgate("display"))

        # swapsubroutine = Qsubroutine(name="swap", gates=[Qgate("swap", qn[2], qn[4])])

        cRysubroutines = []
        for i in range(4):
            cRysubroutines += [cRy(qubitnamea=qn[i+1], qubitnameb=qn[0], n=3-i, r=r)]

        qftsubroutine = QFT(n=4, qubitnames=qn[1:5])

        unexpAsubroutines = []
        for i in reversed(range(4)):
            unexpAsubroutines += [expA(qubitnamea=qn[i+1], qubitnameb=qn[5], qubitnamec=qn[6], sign=-1, n=i)]

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
        subroutines += [reversesubroutine]
        # subroutines += [swapsubroutine]
        subroutines += cRysubroutines
        # subroutines += [swapsubroutine]
        subroutines += [reversesubroutine]
        subroutines += [qftsubroutine]
        subroutines += unexpAsubroutines
        subroutines += [uninitsubroutine]
        subroutines += [resultsubroutine]

        super().__init__(name="Cao2012 Experiment", qubits=qubits, subroutines=subroutines)


class test_expa(Qfunction):

    def __init__(self, m=0, n=0, dorotation=True, noglobalrotation=False):

        name = "exp(2*pi*i*A/16) test\n" \
               "# \n" \
               "# The matrix A is: (1/4)*[15,9,5,-3;9,15,3,-5;5,3,15,-9;-3,-5,-9,15];\n" \
               "# \n" \
               "# The eigenvectors are: [-1;1;1;1], [1;-1;1;1], [1;1;-1;1], [1;1;1;-1] with eigenvalues 1, 2, 4, 8\n" \
               "# \n" \
               "# The value of the first item in one of the vectors is represented by the amplitude of the |00> state, \n" \
               "# The value of the second item by the amplitude of the |01> state, \n" \
               "# The value of the third item by the amplitude of the |10> state, \n" \
               "# The value of the fourth item by the amplitude of the |11> state.\n" \
               "# \n# The eigenvectors can be built with q1 and q2 using:\n" \
               "#    h q1\n" \
               "#    h q2\n" \
               "#    (x q1)\n" \
               "#    [x q2]\n" \
               "#    cz q1,q2\n" \
               "#    [x q2]\n" \
               "#    (x q1)\n" \
               "# The commands in brackets are optional:\n" \
               "# - if none of the bracketed commands are performed, we get [1;1;1;-1], and a rotation of pi/8\n" \
               "# - If only the optional commands in the round brackets are performed, we get [1;-1;1;1], and a rotation of pi/4\n" \
               "# - If only the optional commands in the square brackets are performed, we get [1;1;-1;1], and a rotation of pi/2\n" \
               "# - If both the optional commands in the round and quare brackets are performed, we get [-1;1;1;1], and a rotation of pi"

        if noglobalrotation:
            qubits = 4
        else:
            qubits = 3

        qn = buildnames(qubits, qubitnames="q")

        initgates = []
        if dorotation:
            initgates += [Qgate("x", qn[0])]
        else:
            initgates += [Qgate("#", "x {}".format(qn[0]))]
        initgates += [Qgate("h", qn[1])]
        initgates += [Qgate("h", qn[2])]
        if m == 0 or m == 1:
            initgates += [Qgate("x", qn[1])]
        else:
            initgates += [Qgate("#", "x {}".format(qn[1]))]
        if m == 0 or m == 2:
            initgates += [Qgate("x", qn[2])]
        else:
            initgates += [Qgate("#", "x {}".format(qn[2]))]
        initgates += [Qgate("cz", qn[1], qn[2])]
        if m == 0 or m == 2:
            initgates += [Qgate("x", qn[2])]
        else:
            initgates += [Qgate("#", "x {}".format(qn[2]))]
        if m == 0 or m == 1:
            initgates += [Qgate("x", qn[1])]
        else:
            initgates += [Qgate("#", "x {}".format(qn[1]))]
        if noglobalrotation:
            initgates += [Qgate("x", qn[3])]
        initgates += [Qgate("display")]
        initsubroutine = Qsubroutine(name="init", gates=initgates)

        if noglobalrotation:
            expasubroutine = expA(qubitnamea=qn[0], qubitnameb=qn[1], qubitnamec=qn[2], qubitnamed=qn[3], sign=1, n=n, noglobalrotation=True)
        else:
            expasubroutine = expA(qubitnamea=qn[0], qubitnameb=qn[1], qubitnamec=qn[2], sign=1, n=n)
        expasubroutine.gates.append(Qgate())
        expasubroutine.gates.append(Qgate("display"))

        if noglobalrotation:
            unexpasubroutine = expA(qubitnamea=qn[0], qubitnameb=qn[1], qubitnamec=qn[2], qubitnamed=qn[3], sign=-1, n=n, noglobalrotation=True)
        else:
            unexpasubroutine = expA(qubitnamea=qn[0], qubitnameb=qn[1], qubitnamec=qn[2], sign=-1, n=n)

        resultgates = [Qgate("display"), Qgate("measure"), Qgate("display")]
        resultsubroutine = Qsubroutine(name="result", gates=resultgates)

        subroutines = []
        subroutines += [initsubroutine]
        subroutines += [expasubroutine]
        subroutines += [unexpasubroutine]
        subroutines += [resultsubroutine]

        super().__init__(name=name, qubits=qubits, subroutines=subroutines)
