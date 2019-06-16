#function: Back and forth Quantum Fourier Transform

qubits 3

.init
  x q0
  x q2
  display

.qft
  h q0
  cr q1,q0
  cr q2,q0
  h q1
  cr q2,q1
  h q2

.display
  display

.iqft
  h q2
  cr q2,q1
  cz q2,q1
  h q1
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