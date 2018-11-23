class Qfunction(object):

    def __init__(self, name='', qubits=0, subroutines=[]):
        self.checksubroutines(subroutines)
        self.name = name
        self.qubits = qubits
        if isinstance(subroutines, list):
            self.subroutines = subroutines
        else:
            self.subroutines = [subroutines]

    def __add__(self, subroutines):
        self.checksubroutines(subroutines)
        self.subroutines += subroutines

    def __str__(self):
        string = "#function: {}\n\nqubits {}".format(self.name, self.qubits)
        for subroutine in self.subroutines:
            string += "\n\n"
            string += str(subroutine)
        return string

    def checksubroutines(self, subroutines):
        if type(subroutines) is type([]):
            for subroutine in subroutines:
                if not isinstance(subroutine, Qsubroutine):
                    raise TypeError("Input must be of type Qsubroutine")
            return True
        elif isinstance(subroutines, Qsubroutine):
            return True
        else:
            raise TypeError("Input must be of type Qsubroutine")

class Qsubroutine(object):

    def __init__(self, name='', gates=[], *args, **kwargs):
        self.checkgates(gates)
        self.name = name
        if isinstance(gates, list):
            self.gates = gates
        else:
            self.gates = [gates]

    def __add__(self, gates):
        self.checkgates(gates)
        self.gates = self.gates + gates

    def __str__(self):
        string = ".{}".format(self.name)
        for gate in self.gates:
            string += "\n"
            string += str(gate)
        return string

    def checkgates(self, gates):
        if isinstance(gates, list):
            for gate in gates:
                if not isinstance(gate, Qgate):
                    raise TypeError("Input must be of type Qgate")
            return True
        elif isinstance(gates, Qgate):
            return True
        else:
            raise TypeError("Input must be of type Qgate")

class Qgate(object):

    def __init__(self, name='', *options):
        self.checkgate(name, options)
        self.name = name
        self.options = options

    def __str__(self):
        string = "  "
        string += self.name
        if len(self.options) > 0 and self.name is not "#":
            string += " "
        for i in range(len(self.options)):
            if i is not 0:
                string += ","
            string += self.options[i]
        return string

    def checkgate(self, name, options):
        if not isinstance(name, str):
            raise TypeError("gate name must be string")
        if name in self.gatelst:
            if len(options) is len(self.gatelst[name]):
                for i in range(len(options)):
                    if type(options[i]) is not type(self.gatelst[name][i]):
                        raise TypeError("Expected other data type in gate input")
                return True
            else:
                raise NameError("Gate input size does not match")
        else:
            raise KeyError("Unknown gate type: {}".format(name))

    gatelst = {
        "map":     ("q","q"),
        "measure": (),
        "display": (),
        "#":       ("str"),
        "":        (),
        "h":       ("q"),
        "x":       ("q"),
        "y":       ("q"),
        "z":       ("q"),
        "rx":      ("q",0.0),
        "ry":      ("q",0.0),
        "rz":      ("q",0.0),
        "ph":      ("q"),
        "s":       ("q"),
        "t":       ("q"),
        "tdag":    ("q"),
        "cnot":    ("q","q"),
        "cx":      ("q","q"),
        "toffoli": ("q","q","q"),
        "swap":    ("q","q"),
        "cphase":  ("q","q"),
        "cz":      ("q","q"),
        "cr":      ("q","q"),
        "prepz":   ("q")
    }

def buildnames(n, qubitnames, defaultname="q"):
    if isinstance(qubitnames, list):
        if len(qubitnames) == n:
            qn = qubitnames
        else:
            raise ValueError("Incorrect number of qubit names")
    elif isinstance(qubitnames, str):
        qn = []
        for i in range(n):
            qn += [qubitnames + str(i)]
    else:
        qn = []
        string = defaultname
        for i in range(n):
            qn += [string + str(i)]

    return qn

if __name__ == "__main__":
    testgatecom = Qgate('#', ['blah'])
    testgatenothing = Qgate('', [])
    testgateh = Qgate('h', ['q0'])
    testgatemeas = Qgate('measure', [])
    testsub = Qsubroutine('test', [testgateh,testgateh,testgateh])
    testsub2 = Qsubroutine('test', [testgatecom,testgateh,testgateh,testgateh,testgatenothing,testgatecom,testgateh,testgateh,testgateh,testgatenothing,testgatecom,testgatemeas])
    testfunc = Qfunction('test', 5, [testsub,testsub2])
    print(testfunc)