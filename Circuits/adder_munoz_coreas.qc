#function: Munoz-Coreas Quantum Adder

qubits 13

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
  map q9,b4
  map q10,a4
  map q11,z
  map q11,c
  x ctrl
  display

.add
  cx a1,b1
  cx a2,b2
  cx a3,b3
  cx a4,b4
  cx a3,a4
  cx a2,a3
  cx a1,a2
  toffoli a0,b0,a1
  toffoli a1,b1,a2
  toffoli a2,b2,a3
  toffoli a3,b3,a4
  toffoli a4,b4,c
  display
  toffoli ctrl,c,z
  display
  toffoli a4,b4,c
  toffoli ctrl,a4,b4
  toffoli a3,b3,a4
  toffoli ctrl,a3,b3
  toffoli a2,b2,a3
  toffoli ctrl,a2,b2
  toffoli a1,b1,a2
  toffoli ctrl,a1,b1
  toffoli a0,b0,a1
  toffoli ctrl,a0,b0
  cx a1,a2
  cx a2,a3
  cx a3,a4
  cx a4,b4
  cx a3,b3
  cx a2,b2
  cx a1,b1

.result
  measure
  display