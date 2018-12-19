# from cQASM import *
import QFT
import AdderQFT
import AdderCuccaro
import MultiplierQFT
import Cao2012_Experiment
import os
import subprocess
import math
import numpy as np
import matplotlib.pyplot as plt


def runQX(filename, qubitnum=1, return_res=False, return_raw=False, show_output=False):
    qx_out = subprocess.run(['qx_simulator_0.1_windows_beta.exe', 'Circuits/'+filename+'.qc'], stdout=subprocess.PIPE)
    qx_out_str = qx_out.stdout.decode('utf-8')
    if show_output:
        print(qx_out_str)
    if return_res:
        res = ""
        for i in range(-51, -51 - 4 * qubitnum, -4):
            res += qx_out_str[i]
        # print(res)
        if return_raw:
            return res, qx_out_str
        else:
            return res
    elif return_raw:
        return qx_out_str


if __name__ == "__main__":
    path = "Circuits/"
    if not os.path.exists(path):
        os.makedirs(path)

    run_qft_test = False
    run_add_qft = False
    run_add_qft_ctrl = False
    run_add_cuc = False
    run_add_cuc_ctrl = False
    run_mul_qft = False
    run_Cao2012 = True
    run_test = False

    inp_a = "10101"
    inp_b = "01010"
    inp_ctrl = "1"

    na = len(inp_a)
    nb = len(inp_b)
    nctrl = len(inp_ctrl)

    n = max(na, nb)
    n_tot = na + nb

    if run_qft_test:
        f = open(path + "qftiqft.qc", "w")
        f.write(str(QFT.QFT_iQFTcircuit(inp=inp_a)))
        f.close()

        res_qft_test = runQX('qftiqft', na, return_res=True)
        outp_qft_test = res_qft_test

        print("\n\nQFT-iQFT test:\n\ninput a    = {}\n\noutput a   = {}".format(
            (n-na)*" " + inp_a,
            (n-na)*" " + outp_qft_test))


    if run_add_qft:
        f = open(path + "adder_qft.qc", "w")
        f.write(str(AdderQFT.ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()
        f = open(path + "subtractor_qft.qc", "w")
        f.write(str(AdderQFT.SUBcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()

        res_add_qft = runQX('adder_qft', n_tot, return_res=True)
        outp_a_plus_b_qft = res_add_qft[nb::]

        res_sub_qft = runQX('subtractor_qft', n_tot, return_res=True)
        outp_a_minus_b_qft = res_sub_qft[nb::]

        print("\n\nAdder QFT:\n\ninput a    = {} = {}\ninput b    = {} = {}\n\noutput a+b = {} = {}\noutput a-b = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            outp_a_plus_b_qft, int(outp_a_plus_b_qft, 2),
            outp_a_minus_b_qft, int(outp_a_minus_b_qft, 2)))

    if run_add_qft_ctrl:
        f = open(path + "adder_qft_ctrl.qc", "w")
        f.write(str(AdderQFT.cADDcircuit(inp_a=inp_a, inp_b=inp_b, inp_c=inp_ctrl)))
        f.close()
        f = open(path + "subtractor_qft_ctrl.qc", "w")
        f.write(str(AdderQFT.cSUBcircuit(inp_a=inp_a, inp_b=inp_b, inp_c=inp_ctrl)))
        f.close()

        res_add_qft_ctrl = runQX('adder_qft_ctrl', n_tot + 2, return_res=True)
        outp_a_plus_b_qft_ctrl = res_add_qft_ctrl[nb+2::]

        res_sub_qft_ctrl = runQX('subtractor_qft_ctrl', n_tot + 2, return_res=True)
        outp_a_minus_b_qft_ctrl = res_sub_qft_ctrl[nb+2::]

        print("\n\nAdder QFT Controlled:\n\ninput a    = {} = {}\ninput b    = {} = {}\ninput ctrl = {}\n\noutput a+b = {} = {}\noutput a-b = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            (n-1)*" " + inp_ctrl,
            outp_a_plus_b_qft_ctrl, int(outp_a_plus_b_qft_ctrl, 2),
            outp_a_minus_b_qft_ctrl, int(outp_a_minus_b_qft_ctrl, 2)))

    if run_add_cuc:
        f = open(path + "adder_cuccaro.qc", "w")
        f.write(str(AdderCuccaro.ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()
        f = open(path + "subtractor_cuccaro.qc", "w")
        f.write(str(AdderCuccaro.SUBcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()

        res_add_cuc = runQX('adder_cuccaro', n_tot + 2, return_res=True)
        outp_a_plus_b_cuc = res_add_cuc[-1::-2]

        res_sub_cuc = runQX('subtractor_cuccaro', n_tot + 2, return_res=True)
        outp_a_minus_b_cuc = res_sub_cuc[-1::-2]

        print("\n\nAdder Cuccaro:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\n\noutput a+b = {} = {}\noutput a-b = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            outp_a_plus_b_cuc, int(outp_a_plus_b_cuc, 2),
            outp_a_minus_b_cuc, int(outp_a_minus_b_cuc, 2)))

    if run_add_cuc_ctrl:
        f = open(path + "adder_cuccaro_ctrl.qc", "w")
        f.write(str(AdderCuccaro.cADDcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl)))
        f.close()
        f = open(path + "subtractor_cuccaro_ctrl.qc", "w")
        f.write(str(AdderCuccaro.cSUBcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl)))
        f.close()

        res_add_cuc_ctrl = runQX('adder_cuccaro_ctrl', n_tot + 3, return_res=True)
        outp_a_plus_b_cuc_ctrl = res_add_cuc_ctrl[-1:0:-2]
        res_sub_cuc_ctrl = runQX('subtractor_cuccaro_ctrl', n_tot + 3, return_res=True)
        outp_a_minus_b_cuc_ctrl = res_sub_cuc_ctrl[-1:0:-2]

        print("\n\nAdder Cuccaro Controlled:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\ninput ctrl = {}\n\noutput a+b = {} = {}\noutput a+b = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            n*" " + inp_ctrl,
            outp_a_plus_b_cuc_ctrl, int(outp_a_plus_b_cuc_ctrl, 2),
            outp_a_minus_b_cuc_ctrl, int(outp_a_minus_b_cuc_ctrl, 2)))

    if run_mul_qft:
        f = open(path + "multiplier_qft.qc", "w")
        f.write(str(MultiplierQFT.MULcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()

        res_add_qft = runQX('multiplier_qft', 2*n_tot + 1, return_res=True)
        outp_a_times_b_qft = res_add_qft[na+nb:-1:]

        print("\n\nMultiplier QFT:\n\ninput a    = {} = {}\ninput b    = {} = {}\n\noutput a*b = {} = {}".format(
            nb*" " + inp_a, int(inp_a, 2),
            na*" " + inp_b, int(inp_b, 2),
            outp_a_times_b_qft, int(outp_a_times_b_qft, 2)))

    if run_Cao2012:
        r = 5

        f = open(path + "Cao2012.qc", "w")
        f.write(str(Cao2012_Experiment.Cao2012Experiment(r=r)))
        f.close()

        res_cao2012, raw_cao2012 = runQX('Cao2012', 7, return_res=True, return_raw=True, show_output=True)

        idx = raw_cao2012.rfind('> +')
        idx = raw_cao2012.rfind('> +', 0, idx)
        outp_raw_cao2012 = raw_cao2012[idx-128*38+4:idx+4]

        A = np.zeros((128, 5))
        state = []
        statebin = []
        areal = []
        aimag = []
        atot = []
        for i in range(128):
            state.append(int(raw_cao2012[idx-i*38-1:idx-i*38-8:-1], 2))
            statebin.append(int(raw_cao2012[idx-i*38-1:idx-i*38-8:-1]))
            areal.append(float(raw_cao2012[idx-i*38-29:idx-i*38-20]))
            aimag.append(float(raw_cao2012[idx-i*38-19:idx-i*38-10]))
            atot.append(np.sqrt(areal[-1]**2 + aimag[-1]**2))
            A[127-i][:] = [state[-1], statebin[-1], areal[-1], aimag[-1], atot[-1]]
        A = A[A[:, 0].argsort()]

        print(A)
        print("")
        print(A[64:68, :])

        plt.bar(state, atot)
        # plt.bar(state, 8 + np.log10(atot))
        plt.show()

        outp_bool_cao2012 = res_cao2012[0]
        outp_vec0_cao2012 = res_cao2012[5]
        outp_vec1_cao2012 = res_cao2012[5]

        print("\n\nCao2012 Experiment:\n\ninput r    = {}\n\noutput phi = {}\noutput v0  = {}\noutput v1  = {}".format(
            r,
            outp_bool_cao2012,
            outp_vec0_cao2012,
            outp_vec1_cao2012))

    if run_test:

        f = open(path + "test_expa.qc", "w")
        f.write(str(Cao2012_Experiment.test_expa(m=0, n=0, dorotation=True, noglobalrotation=True)))
        f.close()

        res_cao2012 = runQX('test_expa', 4, show_output=True)
