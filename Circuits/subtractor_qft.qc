#function: QFT Quantum Subtractor

qubits 8

.init
  map q0,b0
  map q1,b1
  map q2,b2
  map q3,b3
  map q4,a0
  map q5,a1
  map q6,a2
  map q7,a3
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

.sub
  cz b0,a0
  cz b1,a1
  swap b1,a1
  cr a1,a0
  cz a1,a0
  swap b1,a1
  cz b2,a2
  swap b2,a2
  cr a2,a1
  cz a2,a1
  cr a2,a0
  swap a2,a1
  cr a1,a0
  swap a2,a1
  cz a2,a0
  swap b2,a2
  cz b3,a3
  swap b3,a3
  cr a3,a2
  cz a3,a2
  cr a3,a1
  swap a3,a2
  cr a2,a1
  swap a3,a2
  cz a3,a1
  cr a3,a0
  swap a3,a2
  cr a2,a0
  swap a3,a2
  swap a3,a1
  cr a1,a0
  swap a3,a1
  cz a3,a0
  swap b3,a3

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