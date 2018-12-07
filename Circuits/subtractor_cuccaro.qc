#function: Quantum Subtractor

qubits 22

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
  map q11,b5
  map q12,a5
  map q13,b6
  map q14,a6
  map q15,b7
  map q16,a7
  map q17,b8
  map q18,a8
  map q19,b9
  map q20,a9
  map q21,z
  x a1
  x a3
  x a5
  x a7
  x a9
  x b0
  x b2
  x b4
  x b6
  x b8
  display

.complement
  x a0
  x a1
  x a2
  x a3
  x a4
  x a5
  x a6
  x a7
  x a8
  x a9

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
  cx a5,b5
  cx a5,a4
  toffoli a4,b5,a5
  cx a6,b6
  cx a6,a5
  toffoli a5,b6,a6
  cx a7,b7
  cx a7,a6
  toffoli a6,b7,a7
  cx a8,b8
  cx a8,a7
  toffoli a7,b8,a8
  cx a9,b9
  cx a9,a8
  toffoli a8,b9,a9
  cx a9,z
  toffoli a8,b9,a9
  cx a9,a8
  cx a8,b9
  toffoli a7,b8,a8
  cx a8,a7
  cx a7,b8
  toffoli a6,b7,a7
  cx a7,a6
  cx a6,b7
  toffoli a5,b6,a6
  cx a6,a5
  cx a5,b6
  toffoli a4,b5,a5
  cx a5,a4
  cx a4,b5
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
  x a5
  x a6
  x a7
  x a8
  x a9
  x b0
  x b1
  x b2
  x b3
  x b4
  x b5
  x b6
  x b7
  x b8
  x b9

.result
  measure
  display