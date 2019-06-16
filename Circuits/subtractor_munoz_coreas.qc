#function: Munoz-Coreas Quantum Subtractor

qubits 11

.init
  map q0,ctrl
  map q1,b0
  map q2,a0
  map q3,b1
  map q4,a1
  map q5,b2
  map q6,a2
  map q7,b3
  map q8,a3
  map q9,z
  map q10,d
  x ctrl
  x a1
  x a3
  x b0
  x b1
  display

.sub
  x a0
  x a1
  x a2
  x a3
  cx a1,b1
  cx a2,b2
  cx a3,b3
  toffoli ctrl,a3,z
  cx a2,a3
  cx a1,a2
  toffoli a0,b0,a1
  toffoli a1,b1,a2
  toffoli a2,b2,a3
  toffoli a3,b3,d
  toffoli ctrl,d,z
  toffoli a3,b3,d
  toffoli ctrl,a3,b3
  toffoli a2,b2,a3
  toffoli ctrl,a2,b2
  toffoli a1,b1,a2
  toffoli ctrl,a1,b1
  toffoli a0,b0,a1
  toffoli ctrl,a0,b0
  cx a1,a2
  cx a2,a3
  cx a3,b3
  cx a2,b2
  cx a1,b1
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