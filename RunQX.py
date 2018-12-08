# from cQASM import *
# from QFT import *
import AdderQFT
import AdderCuccaro
import MultiplierQFT
import os
import subprocess


def runQX(filename, qubitnum):
    qx_out = subprocess.run(['qx_simulator_0.1_windows_beta.exe', 'Circuits/'+filename+'.qc'], stdout=subprocess.PIPE)
    qx_out_str = qx_out.stdout.decode('utf-8')
    res = ""
    for i in range(-51, -51 - 4 * qubitnum, -4):
        res += qx_out_str[i]
    # print(qx_out_str)
    # print(res)
    return res


if __name__ == "__main__":
    path = "Circuits/"
    if not os.path.exists(path):
        os.makedirs(path)

    # f = open(path + "qftiqft.qc", "w")
    # f.write(str(QFT_iQFTcircuit(inp="10101010")))
    # f.close()
    # result = subprocess.run(['qx_simulator_0.1_windows_beta.exe', 'Circuits/qftiqft.qc'], stdout=subprocess.PIPE)
    # outputstr = result.stdout.decode('utf-8')
    # print(outputstr)

    inp_a = "10101"
    inp_b = "01010"
    inp_ctrl = "1"

    run_add_qft = False
    run_add_qft_ctrl = False
    run_add_cuc = False
    run_add_cuc_ctrl = True
    run_mul_qft = False

    na = len(inp_a)
    nb = len(inp_b)
    nctrl = len(inp_ctrl)

    n = max(na, nb)
    n_tot = na + nb

    if run_add_qft:
        f = open(path + "adder_qft.qc", "w")
        f.write(str(AdderQFT.ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()
        f = open(path + "subtractor_qft.qc", "w")
        f.write(str(AdderQFT.SUBcircuit(inp_a=inp_a, inp_b=inp_b)))
        f.close()

        res_add_qft = runQX('adder_qft', n_tot)
        outp_a_plus_b_qft = res_add_qft[nb::]

        res_sub_qft = runQX('subtractor_qft', n_tot)
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

        res_add_qft_ctrl = runQX('adder_qft_ctrl', n_tot + 2)
        outp_a_plus_b_qft_ctrl = res_add_qft_ctrl[nb+2::]

        res_sub_qft_ctrl = runQX('subtractor_qft_ctrl', n_tot + 2)
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

        res_add_cuc = runQX('adder_cuccaro', n_tot + 2)
        outp_a_plus_b_cuc = res_add_cuc[-1::-2]

        res_sub_cuc = runQX('subtractor_cuccaro', n_tot + 2)
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

        res_add_cuc_ctrl = runQX('adder_cuccaro_ctrl', n_tot + 3)
        outp_a_plus_b_cuc_ctrl = res_add_cuc_ctrl[-1:0:-2]
        res_sub_cuc_ctrl = runQX('subtractor_cuccaro_ctrl', n_tot + 3)
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

        res_add_qft = runQX('multiplier_qft', 2*n_tot + 1)
        outp_a_times_b_qft = res_add_qft[na+nb:-1:]

        print("\n\nMultiplier QFT:\n\ninput a    = {} = {}\ninput b    = {} = {}\n\noutput a*b = {} = {}".format(
            nb*" " + inp_a, int(inp_a, 2),
            na*" " + inp_b, int(inp_b, 2),
            outp_a_times_b_qft, int(outp_a_times_b_qft, 2)))
