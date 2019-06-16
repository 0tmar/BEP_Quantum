#function: Controlled QFT Quantum Adder

qubits 10

.init
  map q0,c
  map q1,b0
  map q2,b1
  map q3,b2
  map q4,b3
  map q5,d
  map q6,a0
  map q7,a1
  map q8,a2
  map q9,a3
  x c
  x a0
  x a2
  x b2
  x b3
  display

.qft
  h a0
  cr a1,a0
  cr a2,a0
  cr a3,a0
  h a1
  cr a2,a1
  cr a3,a1
  h a2
  cr a3,a2
  h a3

.cadd
  toffoli c,b0,d
  cz d,a0
  toffoli c,b0,d
  toffoli c,b1,d
  cz d,a1
  swap d,a1
  cr a1,a0
  swap d,a1
  toffoli c,b1,d
  toffoli c,b2,d
  cz d,a2
  swap d,a2
  cr a2,a1
  cr a2,a0
  swap d,a2
  toffoli c,b2,d
  toffoli c,b3,d
  cz d,a3
  swap d,a3
  cr a3,a2
  cr a3,a1
  cr a3,a0
  swap d,a3
  toffoli c,b3,d

.iqft
  h a3
  cr a3,a2
  cz a3,a2
  h a2
  cr a3,a1
  swap a3,a2
  cr a2,a1
  swap a3,a2
  cz a3,a1
  cr a2,a1
  cz a2,a1
  h a1
  cr a3,a0
  swap a3,a2
  cr a2,a0
  swap a3,a2
  swap a3,a1
  cr a1,a0
  swap a3,a1
  cz a3,a0
  cr a2,a0
  swap a2,a1
  cr a1,a0
  swap a2,a1
  cz a2,a0
  cr a1,a0
  cz a1,a0
  h a0

.result
  measure
  display