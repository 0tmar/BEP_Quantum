#function: Controlled Cuccaro Quantum Subtractor

qubits 11

.init
  map q0,ctrl
  map q1,c
  map q2,b0
  map q3,a0
  map q4,b1
  map q5,a1
  map q6,b2
  map q7,a2
  map q8,b3
  map q9,a3
  map q10,z
  x ctrl
  x a1
  x a3
  x b0
  x b1
  display

.csub
  x a0
  x a1
  x a2
  x a3
  toffoli ctrl,a0,b0
  cx a0,c
  toffoli c,b0,a0
  toffoli ctrl,a1,b1
  cx a1,a0
  toffoli a0,b1,a1
  toffoli ctrl,a2,b2
  cx a2,a1
  toffoli a1,b2,a2
  toffoli ctrl,a3,b3
  cx a3,a2
  toffoli a2,b3,a3
  toffoli ctrl,a3,z
  toffoli a2,b3,a3
  cx a3,a2
  toffoli ctrl,a2,b3
  toffoli a1,b2,a2
  cx a2,a1
  toffoli ctrl,a1,b2
  toffoli a0,b1,a1
  cx a1,a0
  toffoli ctrl,a0,b1
  toffoli c,b0,a0
  cx a0,c
  toffoli ctrl,c,b0
  x a0
  x a1
  x a2
  x a3
  cx ctrl,b0
  cx ctrl,b1
  cx ctrl,b2
  cx ctrl,b3

.result
  measure
  display