# from cQASM import *
# from QFT import *
import AdderQFT
import AdderCuccaro
import os
import subprocess


def runQX(filename, qubitnum):
    qx_out = subprocess.run(['qx_simulator_0.1_windows_beta.exe', 'Circuits/'+filename+'.qc'], stdout=subprocess.PIPE)
    qx_out_str = qx_out.stdout.decode('utf-8')
    res = ""
    for i in range(-51, -51 - 4 * qubitnum, -4):
        res += qx_out_str[i]
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

    inp_a = "1010101010"
    inp_b = "0101010101"

    na = len(inp_a)
    nb = len(inp_b)
    n = na + nb

    f = open(path + "adder_qft.qc", "w")
    f.write(str(AdderQFT.ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
    f.close()
    f = open(path + "subtractor_qft.qc", "w")
    f.write(str(AdderQFT.SUBcircuit(inp_a=inp_a, inp_b=inp_b)))
    f.close()
    f = open(path + "adder_cuccaro.qc", "w")
    f.write(str(AdderCuccaro.ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
    f.close()
    f = open(path + "subtractor_cuccaro.qc", "w")
    f.write(str(AdderCuccaro.SUBcircuit(inp_a=inp_a, inp_b=inp_b)))
    f.close()

    res_add_qft = runQX('adder_qft', n)
    outp_a_plus_b_qft = res_add_qft[nb:]

    res_sub_qft = runQX('subtractor_qft', n)
    outp_a_minus_b_qft = res_sub_qft[nb:]

    print("\n\nQFT:\n\ninput a    = {} = {}\ninput b    = {} = {}\n\noutput a+b = {} = {}\noutput a-b = {} = {}".format(
        inp_a, int(inp_a, 2),
        inp_b, int(inp_b, 2),
        outp_a_plus_b_qft, int(outp_a_plus_b_qft, 2),
        outp_a_minus_b_qft, int(outp_a_minus_b_qft, 2)))

    res_add_cuc = runQX('adder_cuccaro', n+2)
    outp_a_plus_b_cuc = res_add_cuc[-1::-2]

    res_sub_cuc = runQX('subtractor_cuccaro', n+2)
    outp_a_minus_b_cuc = res_sub_cuc[-1::-2]

    print("\n\nCuccaro:\n\ninput a    =  {} = {}\ninput b    =  {} = {}\n\noutput a+b = {} = {}\noutput a-b = {} = {}".format(
        inp_a, int(inp_a, 2),
        inp_b, int(inp_b, 2),
        outp_a_plus_b_cuc, int(outp_a_plus_b_cuc, 2),
        outp_a_minus_b_cuc, int(outp_a_minus_b_cuc, 2)))
