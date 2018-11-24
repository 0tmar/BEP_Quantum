from cQASM import *
from QFT import *
from Adder import *
import os
import subprocess

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

    inp_a = "10111"
    inp_b = "10101"
    f = open(path + "adder.qc", "w")
    f.write(str(ADDcircuit(inp_a=inp_a, inp_b=inp_b)))
    f.close()
    f = open(path + "subtractor.qc", "w")
    f.write(str(SUBcircuit(inp_a=inp_a, inp_b=inp_b)))
    f.close()

    addresult = subprocess.run(['qx_simulator_0.1_windows_beta.exe', 'Circuits/adder.qc'], stdout=subprocess.PIPE)
    addoutputstr = addresult.stdout.decode('utf-8')
    #print(addoutputstr)
    addres = ""
    for i in range(-51, -51 - 4 * len(inp_a + inp_b), -4):
        addres += addoutputstr[i]
    outp_a_plus_b = addres[len(inp_b):]

    subresult = subprocess.run(['qx_simulator_0.1_windows_beta.exe', 'Circuits/subtractor.qc'], stdout=subprocess.PIPE)
    suboutputstr = subresult.stdout.decode('utf-8')
    #print(suboutputstr)
    subres = ""
    for i in range(-51, -51 - 4 * len(inp_a + inp_b), -4):
        subres += suboutputstr[i]
    outp_a_minus_b = subres[len(inp_b):]

    print("\ninput a    = {}\ninput b    = {}\n\noutput a+b = {}\noutput a-b = {}".format(inp_a, inp_b, outp_a_plus_b, outp_a_minus_b))
