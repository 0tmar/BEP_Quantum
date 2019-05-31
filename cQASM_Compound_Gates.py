from cQASM import *


def ccZ(qna, qnb, qnc, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' ccZ({},{},{})'.format(qna, qnb, qnc))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 3:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (3), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' ccZ({},{},{})'.format(dcn[0], dcn[1], dcn[2]))]
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
    return gates


def cRx(qna, qnb, theta, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cRx({},{},{})'.format(qna, qnb, theta))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cRx({},{},{})'.format(dcn[0], dcn[1], theta))]
    gates += [Qgate('rx', qnb, theta/4)]
    gates += [Qgate('cz', qna, qnb)]
    gates += [Qgate('rx', qnb, -theta/2)]
    gates += [Qgate('cz', qna, qnb)]
    gates += [Qgate('rx', qnb, theta/4)]
    return gates


def cRy(qna, qnb, theta, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cRx({},{},{})'.format(qna, qnb, theta))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cRx({},{},{})'.format(dcn[0], dcn[1], theta))]
    gates += [Qgate('ry', qnb, theta/4)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('ry', qnb, -theta/2)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('ry', qnb, theta/4)]
    return gates


def cRz(qna, qnb, theta, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cRx({},{},{})'.format(qna, qnb, theta))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cRx({},{},{})'.format(dcn[0], dcn[1], theta))]
    gates += [Qgate('rz', qnb, theta/4)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('rz', qnb, -theta/2)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('rz', qnb, theta/4)]
    return gates


def ccRx(qna, qnb, qnc, theta, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' ccRx({},{},{},{})'.format(qna, qnb, qnc, theta))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 3:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (3), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' ccRx({},{},{},{})'.format(dcn[0], dcn[1], dcn[2], theta))]
    gates += [Qgate("rx", qnc, theta/4)]
    gates += [Qgate("cz", qnb, qnc)]
    gates += [Qgate("rx", qnc, -theta/4)]
    gates += [Qgate("cz", qna, qnb)]
    gates += [Qgate("cz", qna, qnc)]
    gates += [Qgate("rx", qnc, theta/4)]
    gates += [Qgate("cz", qnb, qnc)]
    gates += [Qgate("cz", qna, qnb)]
    gates += [Qgate("cz", qna, qnc)]
    gates += [Qgate("rx", qnc, -theta/4)]
    gates += [Qgate("cz", qna, qnc)]
    return gates


def ccRy(qna, qnb, qnc, theta, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' ccRy({},{},{},{})'.format(qna, qnb, qnc, theta))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 3:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (3), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' ccRy({},{},{},{})'.format(dcn[0], dcn[1], dcn[2], theta))]
    gates += [Qgate("ry", qnc, theta/4)]
    gates += [Qgate("cx", qnb, qnc)]
    gates += [Qgate("ry", qnc, -theta/4)]
    gates += [Qgate("cx", qna, qnb)]
    gates += [Qgate("cx", qna, qnc)]
    gates += [Qgate("ry", qnc, theta/4)]
    gates += [Qgate("cx", qnb, qnc)]
    gates += [Qgate("cx", qna, qnb)]
    gates += [Qgate("cx", qna, qnc)]
    gates += [Qgate("ry", qnc, -theta/4)]
    gates += [Qgate("cx", qna, qnc)]
    return gates


def ccRz(qna, qnb, qnc, theta, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' ccRz({},{},{},{})'.format(qna, qnb, qnc, theta))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 3:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (3), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' ccRz({},{},{},{})'.format(dcn[0], dcn[1], dcn[2], theta))]
    gates += [Qgate("rz", qnc, theta/4)]
    gates += [Qgate("cx", qnb, qnc)]
    gates += [Qgate("rz", qnc, -theta/4)]
    gates += [Qgate("cx", qna, qnb)]
    gates += [Qgate("cx", qna, qnc)]
    gates += [Qgate("rz", qnc, theta/4)]
    gates += [Qgate("cx", qnb, qnc)]
    gates += [Qgate("cx", qna, qnb)]
    gates += [Qgate("cx", qna, qnc)]
    gates += [Qgate("rz", qnc, -theta/4)]
    gates += [Qgate("cx", qna, qnc)]
    return gates


def cH(qna, qnb, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cH({},{})'.format(qna, qnb))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cH({},{})'.format(dcn[0], dcn[1]))]
    gates += [Qgate('t', qna)]
    gates += [Qgate('x', qna)]
    gates += [Qgate('tdag', qna)]
    gates += [Qgate('x', qna)]
    gates += [Qgate('h', qnb)]
    gates += [Qgate('tdag', qnb)]
    gates += [Qgate('tdag', qnb)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('h', qnb)]
    gates += [Qgate('t', qnb)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('t', qnb)]
    gates += [Qgate('h', qnb)]
    gates += [Qgate('s', qnb)]
    gates += [Qgate('x', qnb)]
    return gates


def cS(qna, qnb, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cS({},{})'.format(qna, qnb))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cS({},{})'.format(dcn[0], dcn[1]))]
    gates += [Qgate('t', qna)]
    gates += [Qgate('t', qnb)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('tdag', qnb)]
    gates += [Qgate('cx', qna, qnb)]
    return gates


def cSdag(qna, qnb, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cS({},{})'.format(qna, qnb))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cS({},{})'.format(dcn[0], dcn[1]))]
    gates += [Qgate('tdag', qna)]
    gates += [Qgate('tdag', qnb)]
    gates += [Qgate('cx', qna, qnb)]
    gates += [Qgate('t', qnb)]
    gates += [Qgate('cx', qna, qnb)]
    return gates


def cV(qna, qnb, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cV({},{})'.format(qna, qnb))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cV({},{})'.format(dcn[0], dcn[1]))]
    gates += cH(qna, qnb, add_comment=False)
    gates += cS(qna, qnb, add_comment=False)
    gates += cH(qna, qnb, add_comment=False)
    return gates


def cVdag(qna, qnb, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' cVdag({},{})'.format(qna, qnb))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 2:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs (2), not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cVdag({},{})'.format(dcn[0], dcn[1]))]
    gates += cH(qna, qnb, add_comment=False)
    gates += cSdag(qna, qnb, add_comment=False)
    gates += cH(qna, qnb, add_comment=False)
    return gates


def ncx(qna, qnb, qnc, invert=False, add_comment=True, different_comment_names=None):
    gates = []
    if add_comment:
        gates += [Qgate()]
        if different_comment_names is None:
            gates += [Qgate('#', ' {}x({},...,{},{}), with ancillae {},...,{}'.format(len(qna)*'c', qna[0], qna[-1], qnb[0], qnc[0], qnc[-1]))]
        else:
            dcn = different_comment_names
            if len(dcn) is not 1000:
                raise IndexError(
                    "different_comment_names must be the se size as the amount of inputs, not {}.".format(len(dcn)))
            gates += [Qgate('#', ' cVdag({},{})'.format(dcn[0], dcn[1]))]
    if len(qna) is not len(qnc)+2:
        raise IndexError("Size of qna must be one larger than qnc")
    elif len(qnb) is not 1:
        raise IndexError("Size of qnb must be exactly 1")
    if invert:
        for i in range(len(qna)):
            gates += [Qgate('x', qna[i])]
    gates += [Qgate('toffoli', qna[0], qna[1], qnc[0])]
    for i in range(len(qnc)-1):
        gates += [Qgate('toffoli', qna[i+2], qnc[i], qnc[i+1])]
    gates += [Qgate('toffoli', qna[-1], qnc[-1], qnb[0])]
    for i in reversed(range(len(qnc)-1)):
        gates += [Qgate('toffoli', qna[i+2], qnc[i], qnc[i+1])]
    gates += [Qgate('toffoli', qna[0], qna[1], qnc[0])]
    if invert:
        for i in reversed(range(len(qna))):
            gates += [Qgate('x', qna[i])]
    return gates
