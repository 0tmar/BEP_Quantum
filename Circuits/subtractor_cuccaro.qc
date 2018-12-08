#function: Cuccaro Quantum Subtractor

qubits 12

.init
  map q0,c
  map q1,b0
  map q2,a0
  map q3,b1
  map q4,a1
  map q5,b2
  map q6,a2
  map q7,b3
  map q8,a3
  map q9,b4
  map q10,a4
  map q11,z
  x a0
  x a2
  x a4
  x b1
  x b3
  display

.complement
  x a0
  x a1
  x a2
  x a3
  x a4

.add
  cx a0,b0
  cx a0,c
  toffoli c,b0,a0
  cx a1,b1
  cx a1,a0
  toffoli a0,b1,a1
  cx a2,b2
  cx a2,a1
  toffoli a1,b2,a2
  cx a3,b3
  cx a3,a2
  toffoli a2,b3,a3
  cx a4,b4
  cx a4,a3
  toffoli a3,b4,a4
  cx a4,z
  toffoli a3,b4,a4
  cx a4,a3
  cx a3,b4
  toffoli a2,b3,a3
  cx a3,a2
  cx a2,b3
  toffoli a1,b2,a2
  cx a2,a1
  cx a1,b2
  toffoli a0,b1,a1
  cx a1,a0
  cx a0,b1
  toffoli c,b0,a0
  cx a0,c
  cx c,b0

.complement
  x a0
  x a1
  x a2
  x a3
  x a4
  x b0
  x b1
  x b2
  x b3
  x b4

.result
  measure
  display