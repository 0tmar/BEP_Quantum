# from cQASM import *
import QFT
import AdderQFT
import AdderCuccaro
import MultiplierQFT
import Cao2012_Experiment
import HLL_Linear_Solver
import Number_Inversion
import DivisionThapliyal
import os
import subprocess
import math
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.nan)


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


def build_output_matrix(outp_raw, idx, qubits, lines, do_sort=True):
    linelength = 31 + qubits
    A = np.zeros((lines, 5))
    state = []
    statebin = []
    areal = []
    aimag = []
    atot = []
    for i in range(lines):
        state.append(int(outp_raw[(idx - (i*linelength) - 1):(idx - (i*linelength) - (qubits+1)):(-1)], 2))
        statebin.append(int(outp_raw[(idx - (i*linelength) - 1):(idx - (i*linelength) - (qubits+1)):(-1)]))
        areal.append(float(outp_raw[(idx - (i*linelength) - (qubits+22)):(idx - (i*linelength) - (qubits+13))]))
        aimag.append(float(outp_raw[(idx - (i*linelength) - (qubits+12)):(idx - (i*linelength) - (qubits+3))]))
        atot.append(np.sqrt(areal[-1] ** 2 + aimag[-1] ** 2))
        A[lines - (i+1)][:] = [state[-1], statebin[-1], areal[-1], aimag[-1], atot[-1]]
    if do_sort:
        A = A[A[:, 0].argsort()]
    return A


def find_next_output_matrix(outp_raw, idx=0, return_extra_info=False, do_sort=False):
    L_intro = 46
    L_line_outro = 5
    L_line_no_qubits = 31

    idx0 = outp_raw.find("--------------[quantum state]--------------", idx)

    if idx0 is not -1:

        idx1 = outp_raw.find("> +", idx0)
        idx2 = outp_raw.find("-------------------------------------------", idx1)

        N_qubits = idx1 - idx0 - L_intro + L_line_outro - L_line_no_qubits
        L_line = L_line_no_qubits + N_qubits
        L_states = idx2 - idx0 - L_intro
        N_lines = L_states // L_line

        A = build_output_matrix(outp_raw=outp_raw, idx=idx2-5, qubits=N_qubits, lines=N_lines, do_sort=do_sort)

        if return_extra_info:
            return A, True, idx0, N_qubits, N_lines
        else:
            return A

    else:

        A = np.zeros((0, 5))

        if return_extra_info:
            return A, False, idx0, 0, 0
        else:
            return A


def find_output_matrices(outp_raw, do_sort=False):

    A_lst = []
    idx = -1

    while True:

        A, found_matrix, idx, nqubits, nlines = find_next_output_matrix(outp_raw=outp_raw, idx=idx+1, return_extra_info=True, do_sort=do_sort)

        if found_matrix:
            A_lst += [A]
        else:
            break

    return A_lst


if __name__ == "__main__":
    path = "Circuits/"
    if not os.path.exists(path):
        os.makedirs(path)

    run_qft_test =                  False
    run_add_qft =                   False
    run_add_qft_ctrl =              False
    run_add_cuc =                   False
    run_add_cuc_ctrl =              False
    run_mul_qft =                   False
    run_expa =                      False
    run_Cao2012 =                   False
    run_HLL_test =                  False
    run_numinv_test =               False
    run_division_thapliyal =        True
    run_division_thapliyal_matrix = False
    run_test =                      False

    inp_a = "11111"
    inp_b = "10001"
    inp_ctrl = "1"

    subtype = 'b-a'
    do_overflow = True

    na = len(inp_a)
    nb = len(inp_b)
    nctrl = len(inp_ctrl)

    n = max(na, nb)
    n_tot = na + nb

    if run_qft_test:
        offset_pre = 1
        offset_post = 2

        f = open(path + "qftiqft.qc", "w")
        f.write(str(QFT.QFT_iQFTcircuit(inp=inp_a, offset_pre=offset_pre, offset_post=offset_post)))
        f.close()

        res_qft_test = runQX('qftiqft', na + offset_pre + offset_post, return_res=True, show_output=True)
        outp_qft_test = res_qft_test[offset_pre:(na + offset_pre)]

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
        f.write(str(AdderCuccaro.SUBcircuit(inp_a=inp_a, inp_b=inp_b, subtype=subtype, do_overflow=do_overflow)))
        f.close()

        res_add_cuc = runQX('adder_cuccaro', n_tot + 2, return_res=True)
        outp_a_plus_b_cuc = res_add_cuc[-1::-2]

        res_sub_cuc = runQX('subtractor_cuccaro', n_tot + 2, return_res=True)
        outp_a_minus_b_cuc = res_sub_cuc[-1::-2]

        print("\n\nAdder Cuccaro:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            outp_a_plus_b_cuc, int(outp_a_plus_b_cuc, 2),
            subtype, outp_a_minus_b_cuc, int(outp_a_minus_b_cuc, 2)))

    if run_add_cuc_ctrl:
        f = open(path + "adder_cuccaro_ctrl.qc", "w")
        f.write(str(AdderCuccaro.cADDcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl)))
        f.close()
        f = open(path + "subtractor_cuccaro_ctrl.qc", "w")
        f.write(str(AdderCuccaro.cSUBcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl, subtype=subtype)))
        f.close()

        res_add_cuc_ctrl = runQX('adder_cuccaro_ctrl', n_tot + 3, return_res=True)
        outp_a_plus_b_cuc_ctrl = res_add_cuc_ctrl[-1:0:-2]
        res_sub_cuc_ctrl = runQX('subtractor_cuccaro_ctrl', n_tot + 3, return_res=True)
        outp_a_minus_b_cuc_ctrl = res_sub_cuc_ctrl[-1:0:-2]

        print("\n\nAdder Cuccaro Controlled:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\ninput ctrl = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            n*" " + inp_ctrl,
            outp_a_plus_b_cuc_ctrl, int(outp_a_plus_b_cuc_ctrl, 2),
            subtype, outp_a_minus_b_cuc_ctrl, int(outp_a_minus_b_cuc_ctrl, 2)))

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

    if run_expa:

        f = open(path + "test_expa.qc", "w")
        f.write(str(Cao2012_Experiment.test_expa(m=0, n=0, dorotation=True, noglobalrotation=True)))
        f.close()

        res_expa = runQX('test_expa', 4, show_output=True)

    if run_Cao2012:
        r = 5  # 2^-r factor in ancilla rotation: higher r is higher precision, but lower a chance of |1>. Default r=5
        m = None  # inputs m-th eigenvector of A instead of the multi state (m=0..3 or None) (either m or n must be None)
        n = None  # inputs the n state iso the multi state (n=0: 00, n=1: 01, n=2: 10, n=3: 11)

        do_plot = True

        f = open(path + "Cao2012.qc", "w")
        f.write(str(Cao2012_Experiment.Cao2012Experiment(r=r, m=m, n=n)))
        f.close()

        res_cao2012, raw_cao2012 = runQX('Cao2012', 7, return_res=True, return_raw=True, show_output=True)

        A_lst = find_output_matrices(outp_raw=raw_cao2012, do_sort=False)
        A = A_lst[-2]
        A2 = A_lst[-3]
        A3 = A_lst[-4]

        A = A[A[:, 0].argsort()]
        A2 = A2[A2[:, 0].argsort()]

        print(A3[:, 0:4])
        print("")
        print(A2[:, 0:4])
        print("")
        print(A)
        print("")
        print(A[0:4, :])
        print("")
        print(A[64:68, :])
        print("")
        print(A[0:4, 2])
        print(A[64:68, 2])

        if do_plot:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.bar(A[:, 0], A[:, 2], width=1, align='edge')
            major_xticks = np.arange(0, 129, 16)
            minor_xticks = np.arange(0, 129, 4)
            major_yticks = np.arange(-1, 1.01, 1)
            minor_yticks = np.arange(-1, 1.01, .1)
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
            ax.set_yticks(major_yticks)
            ax.set_yticks(minor_yticks, minor=True)
            ax.grid(which='minor', alpha=0.2)
            ax.grid(which='major', alpha=0.5)
            ax.set_xlim(0, 129)
            ax.set_ylim(-1, 1)
            fig.show()

        outp_bool_cao2012 = res_cao2012[0]
        outp_vec0_cao2012 = res_cao2012[5]
        outp_vec1_cao2012 = res_cao2012[5]

        print("\n\nCao2012 Experiment:\n\ninput r    = {}\n\noutput phi = {}\noutput v0  = {}\noutput v1  = {}".format(
            r,
            outp_bool_cao2012,
            outp_vec0_cao2012,
            outp_vec1_cao2012))

    if run_HLL_test:

        f = open(path + "test_hll.qc", "w")
        f.write(str(HLL_Linear_Solver.EigenvalueInversion_Circuit(n=3, x=1, remove_global_shift=True, save_value=True)))
        f.close()

        raw_HLL = runQX('test_hll', 16, show_output=True, return_raw=True)

        A_lst = find_output_matrices(outp_raw=raw_HLL, do_sort=False)
        A = A_lst[-1]
        print(A[:, 0:4])

    if run_numinv_test:

        f = open(path + "test_numinv.qc", "w")
        f.write(str(Number_Inversion.NUMINVcircuit(inp=inp_a, order=1)))
        f.close()

        raw_HLL = runQX('test_numinv', 3*na+1, show_output=True, return_raw=True)

    if run_division_thapliyal:

        f = open(path + "division_thapliyal.qc", "w")
        f.write(str(DivisionThapliyal.DIVcircuit(inp_n=inp_a, inp_d=inp_b)))
        f.close()

        res_division_thapliyal = runQX('division_thapliyal', 3*n + 1, return_res=True)
        outp_q_division_thapliyal = res_division_thapliyal[(2*n)-1:n-1:-1]
        outp_r_division_thapliyal = res_division_thapliyal[n-1::-1]

        print("\n\nDivision Thapliyal:\n\ninput n        = {} = {}\ninput d        = {} = {}\n\noutput q = n/d = {} = {}\noutput r = n%d = {} = {}\n\ntest q*d + r   = {}*{} + {} = {}\ncorrectness:   {}".format(
            inp_a, int(inp_a, 2),
            inp_b, int(inp_b, 2),
            outp_q_division_thapliyal, int(outp_q_division_thapliyal, 2),
            outp_r_division_thapliyal, int(outp_r_division_thapliyal, 2),
            int(outp_q_division_thapliyal, 2), int(inp_b, 2), int(outp_r_division_thapliyal, 2), int(outp_q_division_thapliyal, 2)*int(inp_b, 2) + int(outp_r_division_thapliyal, 2),
            int(outp_q_division_thapliyal, 2)*int(inp_b, 2) + int(outp_r_division_thapliyal, 2) == int(inp_a, 2)))

    if run_division_thapliyal_matrix:

        N = 5
        n = math.ceil(math.log(N)/math.log(2))
        corr_arr = np.zeros((2**N, 2**N), dtype=bool)

        for i in range(2**N):
            for j in range(2**N):


                f = open(path + "division_thapliyal.qc", "w")
                f.write(str(DivisionThapliyal.DIVcircuit(inp_n=f"{i:0{N}b}", inp_d=f"{j:0{N}b}")))  # Efficiently converts i and j to binary
                f.close()

                res_division_thapliyal = runQX('division_thapliyal', 3*N + 1, return_res=True)
                outp_q = int(res_division_thapliyal[(2*N)-1:N-1:-1], 2)
                outp_r = int(res_division_thapliyal[N-1::-1], 2)

                if not (i==0 or j==0):
                    is_corr = (outp_q == (i//j)) & (outp_r == (i%j))
                else:
                    is_corr = (outp_q == i) & (outp_r == 0)

                if not is_corr:
                    corr_arr[i, j] = is_corr
                    print(f"n={i:>{n}}, d={j:>{n}}:   q={outp_q:>{n}}, r={outp_r:>{n}}")

    if run_test:
        res_test = runQX('test_ccrz', show_output=True)
