# from cQASM import *
import QFT
import AdderQFT
import AdderCuccaro
import AdderMunozCoreas
import MultiplierQFT
import Cao2012Experiment
import NumberInversionCao
import NumberInversionNewtonRaphson
import DivisionThapliyal
import AncillaRotation
import CompleteQLSA
import os
import subprocess
import math
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.nan)


def runQX(filename, qubitnum=1, return_res=False, return_raw=False, show_output=False):
    """Runs the file in 'filename' in the QX simulator
        int  qubitnum:    the amount of qubits that the script uses. Incorrect input yields incorrect answers as output,
        bool return_res:  whether the measured output should be returned as a bitstring,
        bool return_raw:  whether the raw output of the QX simulator should be returned,
        bool show_output: whether the raw output of the QX simulator should be shown."""

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
    """Builds a matrix of the complex state of the qubits in the QX simulator of a 'display' command. This is a 'dumb'
        function, so it relies on another function, 'find_next_output_matrix', to do the searching.
    Outputs:
        np.array A:    matrix containing the information on all found states. The different states are put in different
                       rows. The columns describe, in order, the following properties of each state:
                       0: decimal value of the state (decimal version of next column),
                       1: binary value of the state, with q0 the LSB,
                       2: real part of the complex amplitude of the state,
                       3: imaginary part of the complex amplitude of the state,
                       4: total amplitude of the state.
    Inputs:
        str  outp_raw: the raw output of the QX simulator in which the 'display' command is used,
        int  idx:      index at which the last state of the display ends (with characteristic symbols '> +'),
        int  qubits:   number of qubits used in the script,
        int  lines:    number of lines describing the complete state that is being converted,
        bool do_sort:  whether the states should be sorted with q0 as LSB (by default it is the MSB)."""

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
    """Finds the results of the first 'display' command used earlier in the QX simulation, given some index.
    Outputs:
        np.array A:   the Nx5 matrix containing the info on all found states. When no states are found, then N=0,
        bool     b:   (optional) whether a display command output has been found,
        int      idx: (optional) index at which the last line of the display output stops (idx=-1 if none is found),
        int      n:   (optional) number of qubits used in the circuit (n=0 if none is found),
        int      N:   (optional) number of lines in the found state (N=0 if none is found).
    Inputs:
        str  outp_raw:          the raw output of the QX simulator in which the 'display' command is used,
        int  idx:               index from which the search for a new state display is started,
        bool return_extra_info: if True, it outputs the four extra output values besides the found matrix A,
        bool do_sort:           whether the output matrix should be sorted with q0 as the LSB instead of as the MSB."""

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
    """Finds all outputs of display commands used in the output of a QX simulation run.
    Outputs:
        list A_lst: list containing all output matrixes (Nx5 np.array's) as defined in function 'build_output_matrix',
    Inputs:
        str  outp_raw: the raw output of the QX simulator in which the 'display' command is used,
        bool do_sort:  whether the output matrix should be sorted with q0 as the LSB instead of as the MSB."""

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

    # Booleans determining which quantum circuits are ran
    run_qft_test =                  False
    run_add_qft =                   False
    run_add_qft_ctrl =              False
    run_add_cuc =                   False
    run_add_cuc_ctrl =              False
    run_add_mun_cor =               False
    run_mul_qft =                   False
    run_division_thapliyal =        False
    run_division_thapliyal_matrix = False
    run_ry_c_x_to_the_p =           False
    run_ancilla_rotation =          False
    run_expa =                      False
    run_Cao2012 =                   True
    run_numinv_cao =                False
    run_numinv_newton =             False
    run_complete_qlsa =             False
    run_test =                      False

    # Binary number inputs used in many of the arithmetic circuits
    inp_a = "01011"
    inp_b = "0101"

    # Extra options for some of the circuits
    inp_ctrl = "1"       # For the controlled adders, what the value of the control qubit is; '1' or '0'
    subtype = 'a-b'      # For the subtracters, which subtraction mode is used; 'a-b' or 'b-a'
    do_overflow = True   # For the adders, whether an overflow qubit is used; True or False
    do_ctrl = True       # For the M-C adder, whether controlled version is used; True or False
    c = 1.1              # The c used in the Ry(c*(x^p)) gate and ancilla rotation subroutine
    p = 7                # The p used in the Ry(c*(x^p)) gate
    k = 4                # Determines the maximum order 1+2k in the arcsin approximation in the ancilla rotation

    # Useful lengths for lining out numbers and determining the amount of qubits.
    na = len(inp_a)
    nb = len(inp_b)
    nctrl = len(inp_ctrl)
    n = max(na, nb)
    n_tot = na + nb

    if run_qft_test:
        '''Runs the Quantum Fourier Transform back and forth to test whether it behaves as expected'''

        # Offsets for testing whether it still works with them
        offset_pre = 0
        offset_post = 0

        # Writing the circuit to a file
        f = open(path + "qftiqft.qc", "w")
        f.write(str(QFT.QFT_iQFTcircuit(inp=inp_a, offset_pre=offset_pre, offset_post=offset_post)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        res_qft_test = runQX('qftiqft', na + offset_pre + offset_post, return_res=True, show_output=True)
        outp_qft_test = res_qft_test[offset_pre:(na + offset_pre)]

        # Showing the results
        print("\n\nQFT-iQFT test:\n\ninput a    = {}\n\noutput a   = {}".format(
            (n-na)*" " + inp_a,
            (n-na)*" " + outp_qft_test))

    if run_add_qft:
        '''Runs the Draper QFT Adder and Subtracter'''

        # Writing the circuits to a file
        f = open(path + "adder_qft.qc", "w")
        f.write(str(AdderQFT.ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()
        f = open(path + "subtractor_qft.qc", "w")
        f.write(str(AdderQFT.SUBcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()

        # Running the circuits in the QX simulator and retrieving the results
        res_add_qft = runQX('adder_qft', n_tot, return_res=True)
        outp_a_plus_b_qft = res_add_qft[nb::]
        res_sub_qft = runQX('subtractor_qft', n_tot, return_res=True)
        outp_a_minus_b_qft = res_sub_qft[nb::]

        # Showing the results
        print("\n\nAdder QFT:\n\ninput a    = {} = {}\ninput b    = {} = {}\n\noutput a+b = {} = {}\noutput a-b = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            outp_a_plus_b_qft, int(outp_a_plus_b_qft, 2),
            outp_a_minus_b_qft, int(outp_a_minus_b_qft, 2)))

    if run_add_qft_ctrl:
        '''Runs the controlled version of the Draper QFT Adder and Subtracter'''

        # Writing the circuits to a file
        f = open(path + "adder_qft_ctrl.qc", "w")
        f.write(str(AdderQFT.cADDcircuit(inp_a=inp_a, inp_b=inp_b, inp_c=inp_ctrl)))
        f.close()
        f = open(path + "subtractor_qft_ctrl.qc", "w")
        f.write(str(AdderQFT.cSUBcircuit(inp_a=inp_a, inp_b=inp_b, inp_c=inp_ctrl)))
        f.close()

        # Running the circuits in the QX simulator and retrieving the results
        res_add_qft_ctrl = runQX('adder_qft_ctrl', n_tot + 2, return_res=True)
        outp_a_plus_b_qft_ctrl = res_add_qft_ctrl[nb+2::]
        res_sub_qft_ctrl = runQX('subtractor_qft_ctrl', n_tot + 2, return_res=True)
        outp_a_minus_b_qft_ctrl = res_sub_qft_ctrl[nb+2::]

        # Showing the results
        print("\n\nAdder QFT Controlled:\n\ninput a    = {} = {}\ninput b    = {} = {}\ninput ctrl = {}\n\noutput a+b = {} = {}\noutput a-b = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            (n-1)*" " + inp_ctrl,
            outp_a_plus_b_qft_ctrl, int(outp_a_plus_b_qft_ctrl, 2),
            outp_a_minus_b_qft_ctrl, int(outp_a_minus_b_qft_ctrl, 2)))

    if run_add_cuc:
        '''Runs the Cuccaro Ripple-Carry Adder and Subtracter'''

        # Writing the circuits to a file
        f = open(path + "adder_cuccaro.qc", "w")
        f.write(str(AdderCuccaro.ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()
        f = open(path + "subtractor_cuccaro.qc", "w")
        f.write(str(AdderCuccaro.SUBcircuit(inp_a=inp_a, inp_b=inp_b, subtype=subtype, do_overflow=do_overflow)))
        f.close()

        # Running the circuits in the QX simulator and retrieving the results
        res_add_cuc = runQX('adder_cuccaro', n_tot + 2, return_res=True)
        outp_a_plus_b_cuc = res_add_cuc[-1::-2]
        res_sub_cuc = runQX('subtractor_cuccaro', n_tot + 2, return_res=True)
        outp_a_minus_b_cuc = res_sub_cuc[-1::-2]

        # Showing the results
        print("\n\nAdder Cuccaro:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            outp_a_plus_b_cuc, int(outp_a_plus_b_cuc, 2),
            subtype, outp_a_minus_b_cuc, int(outp_a_minus_b_cuc, 2)))

    if run_add_cuc_ctrl:
        '''Runs the controlled version of the Cuccaro Ripple-Carry Adder and Subtracter'''

        # Writing the circuits to a file
        f = open(path + "adder_cuccaro_ctrl.qc", "w")
        f.write(str(AdderCuccaro.cADDcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl)))
        f.close()
        f = open(path + "subtractor_cuccaro_ctrl.qc", "w")
        f.write(str(AdderCuccaro.cSUBcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl, subtype=subtype)))
        f.close()

        # Running the circuits in the QX simulator and retrieving the results
        res_add_cuc_ctrl = runQX('adder_cuccaro_ctrl', n_tot + 3, return_res=True)
        outp_a_plus_b_cuc_ctrl = res_add_cuc_ctrl[-1:0:-2]
        res_sub_cuc_ctrl = runQX('subtractor_cuccaro_ctrl', n_tot + 3, return_res=True)
        outp_a_minus_b_cuc_ctrl = res_sub_cuc_ctrl[-1:0:-2]

        # Showing the results
        print("\n\nAdder Cuccaro Controlled:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\ninput ctrl = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
            (n-na)*" " + inp_a, int(inp_a, 2),
            (n-nb)*" " + inp_b, int(inp_b, 2),
            n*" " + inp_ctrl,
            outp_a_plus_b_cuc_ctrl, int(outp_a_plus_b_cuc_ctrl, 2),
            subtype, outp_a_minus_b_cuc_ctrl, int(outp_a_minus_b_cuc_ctrl, 2)))

    if run_add_mun_cor:
        '''Runs the Munoz-Coreas Ripple-Carry Adder and Subtracter'''

        # Determines whether an ancilla qubit will be used, and how the output should be read accordingly
        do_ancilla = (do_overflow and do_ctrl)
        if do_overflow and do_ctrl:
            read_start = -2
            read_end = 0
        elif do_ctrl:
            read_start = -2
            read_end = 0
        elif do_overflow:
            read_start = -1
            read_end = None
        else:
            read_start = -2
            read_end = None

        # Writing the circuits to a file
        f = open(path + "adder_munoz_coreas.qc", "w")
        f.write(str(AdderMunozCoreas.ADDcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl, do_overflow=do_overflow, do_ctrl=do_ctrl)))
        f.close()
        f = open(path + "subtractor_munoz_coreas.qc", "w")
        f.write(str(AdderMunozCoreas.SUBcircuit(inp_a=inp_a, inp_b=inp_b, inp_ctrl=inp_ctrl, do_overflow=do_overflow, do_ctrl=do_ctrl, subtype=subtype)))
        f.close()

        # Running the circuits in the QX simulator and retrieving the results
        res_add_mun_cor = runQX('adder_munoz_coreas', n_tot + 1*do_overflow + 1*do_ctrl + 1*do_ancilla, return_res=True, show_output=False)
        outp_a_plus_b_mun_cor = res_add_mun_cor[read_start:read_end:-2]
        res_sub_mun_cor = runQX('subtractor_munoz_coreas', n_tot + 1*do_overflow + 1*do_ctrl + 1*do_ancilla, return_res=True, show_output=False)
        outp_a_minus_b_mun_cor = res_sub_mun_cor[read_start:read_end:-2]

        # Showing the results
        if do_overflow and do_ctrl:
            print("\n\nAdder Munoz-Coreas:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\ninput ctrl = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
                (n-na)*" " + inp_a, int(inp_a, 2),
                (n-nb)*" " + inp_b, int(inp_b, 2),
                n*" " + inp_ctrl,
                outp_a_plus_b_mun_cor, int(outp_a_plus_b_mun_cor, 2),
                subtype, outp_a_minus_b_mun_cor, int(outp_a_minus_b_mun_cor, 2)))
        elif do_overflow:
            print("\n\nAdder Munoz-Coreas:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
                (n-na)*" " + inp_a, int(inp_a, 2),
                (n-nb)*" " + inp_b, int(inp_b, 2),
                outp_a_plus_b_mun_cor, int(outp_a_plus_b_mun_cor, 2),
                subtype, outp_a_minus_b_mun_cor, int(outp_a_minus_b_mun_cor, 2)))
        elif do_ctrl:
            print("\n\nAdder Munoz-Coreas:\n\ninput a    = {} = {}\ninput b    = {} = {}\ninput ctrl = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
                (n-na)*" " + inp_a, int(inp_a, 2),
                (n-nb)*" " + inp_b, int(inp_b, 2),
                (n-1)*" " + inp_ctrl,
                outp_a_plus_b_mun_cor, int(outp_a_plus_b_mun_cor, 2),
                subtype, outp_a_minus_b_mun_cor, int(outp_a_minus_b_mun_cor, 2)))
        else:
            print("\n\nAdder Munoz-Coreas:\n\ninput a    = {} = {}\ninput b    = {} = {}\n\noutput a+b = {} = {}\noutput {} = {} = {}".format(
                (n-na)*" " + inp_a, int(inp_a, 2),
                (n-nb)*" " + inp_b, int(inp_b, 2),
                outp_a_plus_b_mun_cor, int(outp_a_plus_b_mun_cor, 2),
                subtype, outp_a_minus_b_mun_cor, int(outp_a_minus_b_mun_cor, 2)))

    if run_mul_qft:
        '''Runs the straightforward multiplier, using the Draper QFT adder'''

        # Writing the circuit to a file
        f = open(path + "multiplier_qft.qc", "w")
        f.write(str(MultiplierQFT.MULcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        res_add_qft = runQX('multiplier_qft', 2*n_tot + 1, return_res=True)
        outp_a_times_b_qft = res_add_qft[na+nb:-1:]

        # Showing the results
        print("\n\nMultiplier QFT:\n\ninput a    = {} = {}\ninput b    = {} = {}\n\noutput a*b = {} = {}".format(
            nb*" " + inp_a, int(inp_a, 2),
            na*" " + inp_b, int(inp_b, 2),
            outp_a_times_b_qft, int(outp_a_times_b_qft, 2)))

    if run_division_thapliyal:
        '''Runs the Thapliyal interger divider'''

        if inp_b[0] == '1' and na > nb:
            inp_b_temp = "0" + inp_b
            nb_temp = nb + 1
        else:
            inp_b_temp = inp_b
            nb_temp = nb

        # Writing the circuit to a file
        f = open(path + "division_thapliyal.qc", "w")
        f.write(str(DivisionThapliyal.DivUnequalCircuit(inp_n=inp_a, inp_d=inp_b_temp)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        res_division_thapliyal = runQX('division_thapliyal', na + 2*nb_temp, return_res=True)
        outp_q_division_thapliyal = res_division_thapliyal[(na+nb_temp)-1:nb_temp-1:-1]
        outp_r_division_thapliyal = res_division_thapliyal[nb_temp-1::-1]

        # Showing the results
        print("\n\nDivision Thapliyal:\n\ninput n        = {} = {}\ninput d        = {} = {}\n\noutput q = n/d = {} = {}\noutput r = n%d = {} = {}\n\ntest q*d + r   = {}*{} + {} = {}\ncorrectness:   {}".format(
            inp_a, int(inp_a, 2),
            inp_b, int(inp_b, 2),
            outp_q_division_thapliyal, int(outp_q_division_thapliyal, 2),
            outp_r_division_thapliyal, int(outp_r_division_thapliyal, 2),
            int(outp_q_division_thapliyal, 2), int(inp_b, 2), int(outp_r_division_thapliyal, 2), int(outp_q_division_thapliyal, 2)*int(inp_b, 2) + int(outp_r_division_thapliyal, 2),
            int(outp_q_division_thapliyal, 2)*int(inp_b, 2) + int(outp_r_division_thapliyal, 2) == int(inp_a, 2)))

    if run_division_thapliyal_matrix:
        '''Runs the Thapliyal integer divider for N qubits, using all possible inputs to test correctness'''

        # Extra variables
        N = 4                                          # Number of qubits in q, r and d
        n = math.ceil(math.log(N)/math.log(2))         # Maximum size of the decimal outputs
        corr_arr = np.zeros((2**N, 2**N), dtype=bool)  # Array which will be filled with all answers

        # Iteration over all possible divisors and dividends
        for i in range(2**N):
            for j in range(2**N):

                # Writing the circuit to a file
                f = open(path + "division_thapliyal.qc", "w")
                f.write(str(DivisionThapliyal.DIVcircuit(inp_n=f"{i:0{N}b}", inp_d=f"{j:0{N}b}")))
                #                   The strings efficiently convert i and j to binary (of length N, padded with zeroes)
                f.close()

                # Running the circuit in the QX simulator and retrieving the results
                res_division_thapliyal = runQX('division_thapliyal', 3*N + 1, return_res=True)
                outp_q = int(res_division_thapliyal[(2*N)-1:N-1:-1], 2)
                outp_r = int(res_division_thapliyal[N-1::-1], 2)

                # Check for answer correctness
                if not (i==0 or j==0):
                    is_corr = (outp_q == (i//j)) & (outp_r == (i%j))
                else:
                    is_corr = (outp_q == i) & (outp_r == 0)

                # Showing the incorrect results
                if not is_corr:
                    corr_arr[i, j] = is_corr
                    print(f"n={i:>{n}}, d={j:>{n}}:   q={outp_q:>{n}}, r={outp_r:>{n}}")
                    # '>' right aligns, 'n' is the number of characters occupied

    if run_ry_c_x_to_the_p:
        '''Runs the Ry(c*(x^(p)))'''

        # Writing the circuit to a file
        f = open(path + "test_r_y_cx_to_the_p.qc", "w")
        f.write(str(AncillaRotation.Ry_c_x_to_the_p_circuit(inp=inp_a, c=c, p=p)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        raw_Ry_c_x_to_the_p = runQX('test_r_y_cx_to_the_p', n + max(0, p - 2) + 1, show_output=True, return_raw=True)
        A_lst = find_output_matrices(outp_raw=raw_Ry_c_x_to_the_p, do_sort=False)
        prob1 = A_lst[1][1, 2]                  # The probability to find the |0> state
        r = math.asin(prob1)                    # The rotation that has been applied (only correct up to r=pi/2)
        x = int(inp_a, 2)/(2.**(len(inp_a)-1))  # The input value for the c register
        y = c * (x ** p)                        # The value that r should be

        # Showing the results
        print(
            "\n\nRy(c*x^p) rotation:\n\ninput c       = {}\ninput x       = {}\ninput p       = {}\n\noutput sin(r) = {}\noutput r      = {}\n\ntest c*(x^p)  = {}\nrel error     = {}{}".format(
                c,
                x,
                p,
                prob1,
                r,
                y,
                abs(1 - abs(r/y)),
                "\n\nWARNING: THE ROTATION IS LARGER THAN pi/2, MEANING THAT THE INPUT AND OUTPUT WILL NOT MATCH!"*(y > (math.pi/2))))

    if run_ancilla_rotation:
        '''Runs the Ancilla Rotation subroutine of the HHL algorithm'''

        # Writing the circuit to a file
        f = open(path + "test_ancilla_rotation.qc", "w")
        f.write(str(AncillaRotation.AncillaRotationCircuit(inp=inp_a, c=c, k=k)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        raw_ancilla_rotation = runQX('test_ancilla_rotation', n + max(0, 2 * min(n, k) - 1) + 1, show_output=True, return_raw=True)
        A_lst = find_output_matrices(outp_raw=raw_ancilla_rotation, do_sort=False)
        prob1 = A_lst[1][1, 2]                  # The probability to find the |0> state
        x = int(inp_a, 2)/(2.**(len(inp_a)-1))  # The input value for register c
        y = c*x                                 # The desired value for prob1
        arcsin_taylor_factors = np.array([AncillaRotation.arcsin_taylor_factor(p) for p in np.arange(k + 1)])
        powers = np.arange(1, 2 * (k + 1), 2)
        arcsin_cx_approx = np.sum(arcsin_taylor_factors*((c*x)**powers))
        p1_expectation = np.sin(arcsin_cx_approx)

        # Showing the results
        print(
            "\n\nAncilla rotation:\n\ninput c               = {}\ninput x               = {}\ninput m               = {}\n\noutput sin(r)         = {}\n\nexpectation sin(r)    = {}\ntest c*x              = {}\n\nrel expectation error = {}\nrel output error      = {}{}".format(
                c,
                x,
                k,
                prob1,
                p1_expectation,
                y,
                abs((prob1 - p1_expectation) / prob1),
                abs((y - prob1) / y),
                "\n\nWARNING: THE DESIRED PROBABILITY IS LARGER THAN 1, MEANING THAT THE DESIRED OUTPUT CAN NEVER MATCH!"*(y > 1)))

    if run_expa:
        '''Performs the exp(i*A*t) gate with A defined as in the Cao2012 paper'''

        high_res = True
        m = 3
        n = 0
        for m in range(4):
            for n in range(4):

                # Writing the circuit to a file
                f = open(path + "test_expa.qc", "w")
                f.write(str(Cao2012Experiment.test_expa(m=m, n=n, dorotation=True, noglobalrotation=True, highres=high_res)))
                f.close()

                # Running the circuit in the QX simulator and retrieving the results
                raw_expa = runQX('test_expa', 4, return_res=False, return_raw=True, show_output=False)

                A_lst = find_output_matrices(outp_raw=raw_expa, do_sort=True)
                A = A_lst[-3]
                vec = np.array([[(-1.)**(m==0)], [(-1.)**(m==1)], [(-1.)**(m==2)], [(-1.)**(m==3)]])

                outp_vec = np.sum(A[:, 2:4]*vec*np.array([1., 1.j]), 1)
                outp_mean = np.mean(outp_vec, 0)
                outp_delta = outp_vec - outp_mean
                oupt_delta_max = np.max(np.abs(outp_delta))
                rot = np.angle(outp_mean, deg=True)
                rot_exp = 22.5*2**(n+m) % 360
                rot_exp = rot_exp if rot_exp <= 180 else rot_exp-360
                rot_abs_err = abs(abs(rot)-abs(rot_exp))
                rot_rel_err = abs(rot_abs_err/rot)

                print("\n\nTest exp(i*A*t_0*2^n)*v_m for Cao matrix:\n\ninput n        = {}\ninput m        = {}\ninput ctrl     = {}\n\noutput delta   = {}\noutput r (deg) = {}\n\nexpected r     = {}\n\nabs error      = {}\nrel error      = {}".format(
                    n,
                    m,
                    1,
                    oupt_delta_max,
                    rot,
                    rot_exp,
                    rot_abs_err,
                    rot_rel_err))

    if run_Cao2012:
        '''Runs the proprietary HHL circuit from Cao2012'''

        r = 5     # 2^-r factor in ancilla rotation: higher r is higher precision, but lower chance of |1>. Default r=5
        k = None  # inputs mth eigenvector of A instead of the multi state (m=0..3 or None) (either m or n must be None)
        n = None  # inputs the n state iso the multi state (n=0: 00, n=1: 01, n=2: 10, n=3: 11)

        do_plot = True
        do_save = False

        # Writing the circuit to a file
        f = open(path + "Cao2012.qc", "w")
        f.write(str(Cao2012Experiment.Cao2012Experiment(r=r, m=k, n=n)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        res_cao2012, raw_cao2012 = runQX('Cao2012', 7, return_res=True, return_raw=True, show_output=True)
        A_lst = find_output_matrices(outp_raw=raw_cao2012, do_sort=False)
        A = A_lst[-2]
        A2 = A_lst[-3]
        A3 = A_lst[-4]

        A = A[A[:, 0].argsort()]
        A2 = A2[A2[:, 0].argsort()]

        # # Printing some inside information
        # print(A3[:, 0:4])  # Complete state right before the iQFT
        # print("")
        # print(A2[:, 0:4])  # Complete state right after the iQFT
        # print("")
        # print(A)  # Complete final state
        # print("")
        # print(A[0:4, :])  # Complete final state minus the total amplitude of each state
        # print("")
        # print(A[64:68, :])  # States of the qubits carrying the information of the output vector
        # print("")
        # print(A[0:4, 2])  # Real amplitudes of the garbage states
        # print(A[64:68, 2])  # Real amplitudes of the qubits carrying the information of the output vector
        # print(7*A[64:68, 2]/A[65, 2])  # Relative version of previous line

        # Plotting the results
        if do_plot:
            fig = plt.figure(1)
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
            if do_save:
                fig.savefig(f'Images/Cao2012_output_r{r:d}.png', bbox_inches='tight')

        # # The (irrelevant) measured outputs of the circuit
        # outp_bool_cao2012 = res_cao2012[0]
        # outp_vec0_cao2012 = res_cao2012[5]
        # outp_vec1_cao2012 = res_cao2012[5]
        #
        # # Showing the (ir)relevant measured qubit outputs
        # print("\n\nCao2012 Experiment:\n\ninput r    = {}\n\noutput phi = {}\noutput v0  = {}\noutput v1  = {}".format(
        #     r,
        #     outp_bool_cao2012,
        #     outp_vec0_cao2012,
        #     outp_vec1_cao2012))

        outp_garbage = A[0:4, 2]
        outp_x = A[64:68, 2]
        outp_x_ideal = np.array([-1, 7, 11, 13])*math.pi/(2**(2+r))
        outp_prob_x = sum(outp_x**2)

        np.set_printoptions(precision=5)    # Number of decimals
        np.set_printoptions(suppress=True)  # Whether to suppress scientific notation
        print("\n\nCao2012 Experiment:\n\ninput r        = {}\n\noutput x       = {}\noutput x ideal = {}\n\noutput x rel   = {}\n\nabs error      = {}\nrel error      = {}".format(
                r,
                outp_x,
                outp_x_ideal,
                -outp_x/outp_x[0],
                abs(outp_x-outp_x_ideal),
                1 - outp_x/outp_x_ideal))

        np.set_printoptions(suppress=True)  # Whether to suppress scientific notation
        print("\nP(|1>) = {}".format(
            outp_prob_x))

    if run_numinv_cao:
        '''Runs the Cao number inversion algorithm'''

        # Writing the circuit to a file
        f = open(path + "test_hll.qc", "w")
        f.write(str(NumberInversionCao.EigenvalueInversion_Circuit(
            n=3,
            x=1,
            remove_global_shift=True,
            save_value=True)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        raw_HLL = runQX('test_hll', 16, show_output=True, return_raw=True)
        A_lst = find_output_matrices(outp_raw=raw_HLL, do_sort=False)
        A = A_lst[-1]
        print(A[:, 0:4])

    if run_numinv_newton:
        '''Performs the number inversion based on Newton-Raphson iteration'''

        # Writing the circuit to a file
        f = open(path + "test_numinv.qc", "w")
        f.write(str(NumberInversionNewtonRaphson.NUMINVcircuit(inp=inp_a, order=1)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        res_numinv_newton_raphson, raw_numinv_newton_raphson = runQX('test_numinv', 3*na+2, show_output=True, return_res=True, return_raw=True)
        A_lst = find_output_matrices(outp_raw=raw_numinv_newton_raphson, do_sort=False)
        outp_x0_bin = int(A_lst[1][0, 1])
        outp_x0_str = str(outp_x0_bin)
        outp_x0_str = (3*na+2 - len(outp_x0_str))*"0" + outp_x0_str
        x0_str = outp_x0_str[na:3*na-1]
        x0 = int(x0_str, 2)/2**(2*na - 2)
        x1_str = res_numinv_newton_raphson[na+1:3*na]
        x1 = int(x1_str, 2)/2**(2*na - 2)
        a = int(inp_a, 2)
        x0_th = 2**(-math.floor(math.log2(a)))
        x0_th_int = int(x0_th*2**(2*na - 2))
        x1_th = 2*x0_th - a*(x0_th**2)
        x1_th_int = int(x1_th*2**(2*na - 2))
        a_inv = 1/a
        a_inv_int = int(2**(2*na - 2)/a)

        print("\n\nInversion Newton-Raphson:\n\ninput a      = {} = {}\n\noutput x0    = {} = {}\noutput x1    = {} = {}\n\nexpected x0  = {} = {}\nexpected x1  = {} = {}\n\noutput ideal = {} = {}".format(
            (na-1) * " " + inp_a, a,
            x0_str, x0,
            x1_str, x1,
            f'{x0_th_int:0{2*na-1}b}', x0_th,
            f'{x1_th_int:0{2*na-1}b}', x1_th,
            f'{a_inv_int:0{2*na-1}b}', a_inv))

    if run_complete_qlsa:
        '''Runs a complete QLSA using the Cao matrix'''

        # Defining extra variables
        n_eig_inv = 4
        m_anc_rot = 9
        r_anc_rot = 5
        input_eigenvector = None
        input_n_state = None
        do_plot = True
        do_save = False

        # Writing the circuit to a file
        f = open(path + "complete_qlsa_cao_matrix.qc", "w")
        f.write(str(CompleteQLSA.CompleteQlsaWithCaoMatrix(
            n_eig_inv=n_eig_inv,
            k_anc_rot=m_anc_rot,
            r_anc_rot=r_anc_rot-2,
            input_eigenvector=input_eigenvector,
            input_n_state=input_n_state)))
        f.close()

        # Running the circuit in the QX simulator and retrieving the results
        raw_qlsa = runQX(
            'complete_qlsa_cao_matrix',
            2 + 4 + (4 + n_eig_inv + 3) + (1 + max(0, min(n_eig_inv+1, 1+2*m_anc_rot) - 2)),
            show_output=True,
            return_raw=True)
        A_lst = find_output_matrices(outp_raw=raw_qlsa, do_sort=False)
        A = A_lst[-2]
        A = A[A[:, 0].argsort()]

        outcome_rel = -A[64:68, 2]/A[64, 2]
        print(outcome_rel)

        # Plotting the results
        if do_plot:
            fig = plt.figure(1)
            ax = fig.add_subplot(1, 1, 1)
            # ax.bar(A[:, 0], A[:, 2], width=1, align='edge')
            ax.bar(np.arange(2**7), A[:, 2], width=1, align='edge')
            major_xticks = np.arange(0, 2**7 + 1, 16)
            minor_xticks = np.arange(0, 2**7 + 1, 4)
            major_yticks = np.arange(-1, 1.01, 1)
            minor_yticks = np.arange(-1, 1.01, .1)
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
            ax.set_yticks(major_yticks)
            ax.set_yticks(minor_yticks, minor=True)
            ax.grid(which='minor', alpha=0.2)
            ax.grid(which='major', alpha=0.5)
            ax.set_xlim(0, 2**7 + 1)
            ax.set_ylim(-1, 1)
            fig.show()
            if do_save:
                fig.savefig(f'Images/Cao2012_output_r{r_anc_rot:d}_m{m_anc_rot:d}.png', bbox_inches='tight')

        outp_garbage = A[0:4, 2]
        outp_x = A[64:68, 2]
        outp_x_ideal = np.array([-1, 7, 11, 13])*math.pi/(2**(2+r_anc_rot))
        outp_prob_x = sum(outp_x**2)

        np.set_printoptions(precision=5)    # Number of decimals
        np.set_printoptions(suppress=True)  # Whether to suppress scientific notation
        print("\n\nCao2012 Experiment general circuit:\n\ninput r        = {}\ninput m        = {}\n\noutput x       = {}\noutput x ideal = {}\n\noutput x rel   = {}\n\nabs error      = {}\nrel error      = {}".format(
            r_anc_rot,
            m_anc_rot,
            outp_x,
            outp_x_ideal,
            -outp_x/outp_x[0],
            abs(outp_x-outp_x_ideal),
            1 - outp_x/outp_x_ideal))

        np.set_printoptions(suppress=True)  # Whether to suppress scientific notation
        print("\nP(|1>) = {}".format(
            outp_prob_x))

    if run_test:
        '''Runs any circuit requested, specifically to test the self-developed controlled gates'''

        # Running the circuit in the QX simulator and retrieving the results
        res_test = runQX('test_ccrz', show_output=True)
