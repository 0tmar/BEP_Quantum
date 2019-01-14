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
