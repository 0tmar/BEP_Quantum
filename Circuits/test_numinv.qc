#function: Cao Number Inversion Estimator

qubits 14

.init
  map q0,a0
  map q1,a1
  map q2,a2
  map q3,a3
  map q4,b0
  map q5,b1
  map q6,b2
  map q7,b3
  map q8,b4
  map q9,b5
  map q10,b6
  map q11,b7
  map q12,c
  map q13,d
  x a1
  x a2
  display

.powerinvert
  cx a0,b3
  cx a0,c
  x c
  toffoli a1,c,b2
  x c
  cx b2,c
  x c
  toffoli a2,c,b1
  x c
  cx b1,c
  x c
  toffoli a3,c,b0
  x c
  cx b0,c
  x a0
  x a1
  x a2
  x a3
  toffoli a0,a1,b4
  toffoli a2,b4,b5
  toffoli a3,b5,c
  toffoli a2,b4,b5
  toffoli a0,a1,b4
  x a3
  x a2
  x a1
  x a0
  x c
  display

.cnot
  cx b3,c

.qft
  h b3
  cr b4,b3
  cr b5,b3
  cr b6,b3
  cr b7,b3
  h b4
  cr b5,b4
  cr b6,b4
  cr b7,b4
  h b5
  cr b6,b5
  cr b7,b5
  h b6
  cr b7,b6
  h b7

.csub
  # cSUB on b3,...,b7 by a0,...,a3 controlled by c
  toffoli c,a0,d
  cz d,b4
  swap d,b4
  cr b4,b3
  cz b4,b3
  swap d,b4
  toffoli c,a0,d
  toffoli c,a1,d
  cz d,b5
  swap d,b5
  cr b5,b4
  cz b5,b4
  cr b5,b3
  swap b5,b4
  cr b4,b3
  swap b5,b4
  cz b5,b3
  swap d,b5
  toffoli c,a1,d
  toffoli c,a2,d
  cz d,b6
  swap d,b6
  cr b6,b5
  cz b6,b5
  cr b6,b4
  swap b6,b5
  cr b5,b4
  swap b6,b5
  cz b6,b4
  cr b6,b3
  swap b6,b5
  cr b5,b3
  swap b6,b5
  swap b6,b4
  cr b4,b3
  swap b6,b4
  cz b6,b3
  swap d,b6
  toffoli c,a2,d
  toffoli c,a3,d
  cz d,b7
  swap d,b7
  cr b7,b6
  cz b7,b6
  cr b7,b5
  swap b7,b6
  cr b6,b5
  swap b7,b6
  cz b7,b5
  cr b7,b4
  swap b7,b6
  cr b6,b4
  swap b7,b6
  swap b7,b5
  cr b5,b4
  swap b7,b5
  cz b7,b4
  cr b7,b3
  swap b7,b6
  cr b6,b3
  swap b7,b6
  swap b7,b5
  cr b5,b3
  swap b7,b5
  swap b7,b4
  cr b4,b3
  swap b7,b4
  cz b7,b3
  swap d,b7
  toffoli c,a3,d

.iqft
  h b7
  cr b7,b6
  cz b7,b6
  h b6
  cr b7,b5
  swap b7,b6
  cr b6,b5
  swap b7,b6
  cz b7,b5
  cr b6,b5
  cz b6,b5
  h b5
  cr b7,b4
  swap b7,b6
  cr b6,b4
  swap b7,b6
  swap b7,b5
  cr b5,b4
  swap b7,b5
  cz b7,b4
  cr b6,b4
  swap b6,b5
  cr b5,b4
  swap b6,b5
  cz b6,b4
  cr b5,b4
  cz b5,b4
  h b4
  cr b7,b3
  swap b7,b6
  cr b6,b3
  swap b7,b6
  swap b7,b5
  cr b5,b3
  swap b7,b5
  swap b7,b4
  cr b4,b3
  swap b7,b4
  cz b7,b3
  cr b6,b3
  swap b6,b5
  cr b5,b3
  swap b6,b5
  swap b6,b4
  cr b4,b3
  swap b6,b4
  cz b6,b3
  cr b5,b3
  swap b5,b4
  cr b4,b3
  swap b5,b4
  cz b5,b3
  cr b4,b3
  cz b4,b3
  h b3
  display
  prepz c

.cnot
  cx b2,c

.qft
  h b2
  cr b3,b2
  cr b4,b2
  cr b5,b2
  h b3
  cr b4,b3
  cr b5,b3
  h b4
  cr b5,b4
  h b5

.csub
  # cSUB on b2,...,b5 by a1,...,a3 controlled by c
  toffoli c,a1,d
  cz d,b3
  swap d,b3
  cr b3,b2
  cz b3,b2
  swap d,b3
  toffoli c,a1,d
  toffoli c,a2,d
  cz d,b4
  swap d,b4
  cr b4,b3
  cz b4,b3
  cr b4,b2
  swap b4,b3
  cr b3,b2
  swap b4,b3
  cz b4,b2
  swap d,b4
  toffoli c,a2,d
  toffoli c,a3,d
  cz d,b5
  swap d,b5
  cr b5,b4
  cz b5,b4
  cr b5,b3
  swap b5,b4
  cr b4,b3
  swap b5,b4
  cz b5,b3
  cr b5,b2
  swap b5,b4
  cr b4,b2
  swap b5,b4
  swap b5,b3
  cr b3,b2
  swap b5,b3
  cz b5,b2
  swap d,b5
  toffoli c,a3,d

.iqft
  h b5
  cr b5,b4
  cz b5,b4
  h b4
  cr b5,b3
  swap b5,b4
  cr b4,b3
  swap b5,b4
  cz b5,b3
  cr b4,b3
  cz b4,b3
  h b3
  cr b5,b2
  swap b5,b4
  cr b4,b2
  swap b5,b4
  swap b5,b3
  cr b3,b2
  swap b5,b3
  cz b5,b2
  cr b4,b2
  swap b4,b3
  cr b3,b2
  swap b4,b3
  cz b4,b2
  cr b3,b2
  cz b3,b2
  h b2
  display
  prepz c

.cnot
  cx b1,c

.qft
  h b1
  cr b2,b1
  cr b3,b1
  h b2
  cr b3,b2
  h b3

.csub
  # cSUB on b1,...,b3 by a2,...,a3 controlled by c
  toffoli c,a2,d
  cz d,b2
  swap d,b2
  cr b2,b1
  cz b2,b1
  swap d,b2
  toffoli c,a2,d
  toffoli c,a3,d
  cz d,b3
  swap d,b3
  cr b3,b2
  cz b3,b2
  cr b3,b1
  swap b3,b2
  cr b2,b1
  swap b3,b2
  cz b3,b1
  swap d,b3
  toffoli c,a3,d

.iqft
  h b3
  cr b3,b2
  cz b3,b2
  h b2
  cr b3,b1
  swap b3,b2
  cr b2,b1
  swap b3,b2
  cz b3,b1
  cr b2,b1
  cz b2,b1
  h b1
  display
  prepz c

.cnot
  cx b0,c

.qft
  h b0
  cr b1,b0
  h b1

.csub
  # cSUB on b0,...,b1 by a3,...,a3 controlled by c
  toffoli c,a3,d
  cz d,b1
  swap d,b1
  cr b1,b0
  cz b1,b0
  swap d,b1
  toffoli c,a3,d

.iqft
  h b1
  cr b1,b0
  cz b1,b0
  h b0
  display
  prepz c

.result
  measure
  display