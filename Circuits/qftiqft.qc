#function: Back and forth Quantum Fourier Transform

qubits 5

.init
  x q0
  x q2
  x q4
  display

.qft
  h q0
  cr q1,q0
  cr q2,q0
  cr q3,q0
  cr q4,q0
  h q1
  cr q2,q1
  cr q3,q1
  cr q4,q1
  h q2
  cr q3,q2
  cr q4,q2
  h q3
  cr q4,q3
  h q4

.iqft
  h q4
  cr q4,q3
  cz q4,q3
  h q3
  cr q4,q2
  swap q4,q3
  cr q3,q2
  swap q4,q3
  cz q4,q2
  cr q3,q2
  cz q3,q2
  h q2
  cr q4,q1
  swap q4,q3
  cr q3,q1
  swap q4,q3
  swap q4,q2
  cr q2,q1
  swap q4,q2
  cz q4,q1
  cr q3,q1
  swap q3,q2
  cr q2,q1
  swap q3,q2
  cz q3,q1
  cr q2,q1
  cz q2,q1
  h q1
  cr q4,q0
  swap q4,q3
  cr q3,q0
  swap q4,q3
  swap q4,q2
  cr q2,q0
  swap q4,q2
  swap q4,q1
  cr q1,q0
  swap q4,q1
  cz q4,q0
  cr q3,q0
  swap q3,q2
  cr q2,q0
  swap q3,q2
  swap q3,q1
  cr q1,q0
  swap q3,q1
  cz q3,q0
  cr q2,q0
  swap q2,q1
  cr q1,q0
  swap q2,q1
  cz q2,q0
  cr q1,q0
  cz q1,q0
  h q0

.result
  measure
  display