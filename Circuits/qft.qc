#function: Quantum Fourier Transform

qubits 8

.init
  x q0
  x q2
  x q4
  x q6
  display

.qft
  h q0
  cr q1,q0
  cr q2,q0
  cr q3,q0
  cr q4,q0
  cr q5,q0
  cr q6,q0
  cr q7,q0
  h q1
  cr q2,q1
  cr q3,q1
  cr q4,q1
  cr q5,q1
  cr q6,q1
  cr q7,q1
  h q2
  cr q3,q2
  cr q4,q2
  cr q5,q2
  cr q6,q2
  cr q7,q2
  h q3
  cr q4,q3
  cr q5,q3
  cr q6,q3
  cr q7,q3
  h q4
  cr q5,q4
  cr q6,q4
  cr q7,q4
  h q5
  cr q6,q5
  cr q7,q5
  h q6
  cr q7,q6
  h q7
  swap q0,q7
  swap q1,q6
  swap q2,q5
  swap q3,q4

.result
  measure
  display