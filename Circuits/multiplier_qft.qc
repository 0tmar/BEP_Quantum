#function: Cuccaro Quantum Multiplier

qubits 17

.init
  map q0,a0
  map q1,a1
  map q2,a2
  map q3,a3
  map q4,b0
  map q5,b1
  map q6,b2
  map q7,b3
  map q8,c0
  map q9,c1
  map q10,c2
  map q11,c3
  map q12,c4
  map q13,c5
  map q14,c6
  map q15,c7
  map q16,z
  x a0
  x a2
  x b2
  x b3
  display

.qft
  h c0
  cr c1,c0
  cr c2,c0
  cr c3,c0
  cr c4,c0
  cr c5,c0
  cr c6,c0
  cr c7,c0
  h c1
  cr c2,c1
  cr c3,c1
  cr c4,c1
  cr c5,c1
  cr c6,c1
  cr c7,c1
  h c2
  cr c3,c2
  cr c4,c2
  cr c5,c2
  cr c6,c2
  cr c7,c2
  h c3
  cr c4,c3
  cr c5,c3
  cr c6,c3
  cr c7,c3
  h c4
  cr c5,c4
  cr c6,c4
  cr c7,c4
  h c5
  cr c6,c5
  cr c7,c5
  h c6
  cr c7,c6
  h c7

.cadd
  toffoli a3,b0,z
  cz z,c4
  swap z,c4
  cr c4,c3
  cr c4,c2
  cr c4,c1
  cr c4,c0
  swap z,c4
  toffoli a3,b0,z
  toffoli a3,b1,z
  cz z,c5
  swap z,c5
  cr c5,c4
  cr c5,c3
  cr c5,c2
  cr c5,c1
  cr c5,c0
  swap z,c5
  toffoli a3,b1,z
  toffoli a3,b2,z
  cz z,c6
  swap z,c6
  cr c6,c5
  cr c6,c4
  cr c6,c3
  cr c6,c2
  cr c6,c1
  cr c6,c0
  swap z,c6
  toffoli a3,b2,z
  toffoli a3,b3,z
  cz z,c7
  swap z,c7
  cr c7,c6
  cr c7,c5
  cr c7,c4
  cr c7,c3
  cr c7,c2
  cr c7,c1
  cr c7,c0
  swap z,c7
  toffoli a3,b3,z

.iqft
  h c7
  cr c7,c6
  cz c7,c6
  h c6
  cr c7,c5
  swap c7,c6
  cr c6,c5
  swap c7,c6
  cz c7,c5
  cr c6,c5
  cz c6,c5
  h c5
  cr c7,c4
  swap c7,c6
  cr c6,c4
  swap c7,c6
  swap c7,c5
  cr c5,c4
  swap c7,c5
  cz c7,c4
  cr c6,c4
  swap c6,c5
  cr c5,c4
  swap c6,c5
  cz c6,c4
  cr c5,c4
  cz c5,c4
  h c4
  cr c7,c3
  swap c7,c6
  cr c6,c3
  swap c7,c6
  swap c7,c5
  cr c5,c3
  swap c7,c5
  swap c7,c4
  cr c4,c3
  swap c7,c4
  cz c7,c3
  cr c6,c3
  swap c6,c5
  cr c5,c3
  swap c6,c5
  swap c6,c4
  cr c4,c3
  swap c6,c4
  cz c6,c3
  cr c5,c3
  swap c5,c4
  cr c4,c3
  swap c5,c4
  cz c5,c3
  cr c4,c3
  cz c4,c3
  h c3
  cr c7,c2
  swap c7,c6
  cr c6,c2
  swap c7,c6
  swap c7,c5
  cr c5,c2
  swap c7,c5
  swap c7,c4
  cr c4,c2
  swap c7,c4
  swap c7,c3
  cr c3,c2
  swap c7,c3
  cz c7,c2
  cr c6,c2
  swap c6,c5
  cr c5,c2
  swap c6,c5
  swap c6,c4
  cr c4,c2
  swap c6,c4
  swap c6,c3
  cr c3,c2
  swap c6,c3
  cz c6,c2
  cr c5,c2
  swap c5,c4
  cr c4,c2
  swap c5,c4
  swap c5,c3
  cr c3,c2
  swap c5,c3
  cz c5,c2
  cr c4,c2
  swap c4,c3
  cr c3,c2
  swap c4,c3
  cz c4,c2
  cr c3,c2
  cz c3,c2
  h c2
  cr c7,c1
  swap c7,c6
  cr c6,c1
  swap c7,c6
  swap c7,c5
  cr c5,c1
  swap c7,c5
  swap c7,c4
  cr c4,c1
  swap c7,c4
  swap c7,c3
  cr c3,c1
  swap c7,c3
  swap c7,c2
  cr c2,c1
  swap c7,c2
  cz c7,c1
  cr c6,c1
  swap c6,c5
  cr c5,c1
  swap c6,c5
  swap c6,c4
  cr c4,c1
  swap c6,c4
  swap c6,c3
  cr c3,c1
  swap c6,c3
  swap c6,c2
  cr c2,c1
  swap c6,c2
  cz c6,c1
  cr c5,c1
  swap c5,c4
  cr c4,c1
  swap c5,c4
  swap c5,c3
  cr c3,c1
  swap c5,c3
  swap c5,c2
  cr c2,c1
  swap c5,c2
  cz c5,c1
  cr c4,c1
  swap c4,c3
  cr c3,c1
  swap c4,c3
  swap c4,c2
  cr c2,c1
  swap c4,c2
  cz c4,c1
  cr c3,c1
  swap c3,c2
  cr c2,c1
  swap c3,c2
  cz c3,c1
  cr c2,c1
  cz c2,c1
  h c1
  cr c7,c0
  swap c7,c6
  cr c6,c0
  swap c7,c6
  swap c7,c5
  cr c5,c0
  swap c7,c5
  swap c7,c4
  cr c4,c0
  swap c7,c4
  swap c7,c3
  cr c3,c0
  swap c7,c3
  swap c7,c2
  cr c2,c0
  swap c7,c2
  swap c7,c1
  cr c1,c0
  swap c7,c1
  cz c7,c0
  cr c6,c0
  swap c6,c5
  cr c5,c0
  swap c6,c5
  swap c6,c4
  cr c4,c0
  swap c6,c4
  swap c6,c3
  cr c3,c0
  swap c6,c3
  swap c6,c2
  cr c2,c0
  swap c6,c2
  swap c6,c1
  cr c1,c0
  swap c6,c1
  cz c6,c0
  cr c5,c0
  swap c5,c4
  cr c4,c0
  swap c5,c4
  swap c5,c3
  cr c3,c0
  swap c5,c3
  swap c5,c2
  cr c2,c0
  swap c5,c2
  swap c5,c1
  cr c1,c0
  swap c5,c1
  cz c5,c0
  cr c4,c0
  swap c4,c3
  cr c3,c0
  swap c4,c3
  swap c4,c2
  cr c2,c0
  swap c4,c2
  swap c4,c1
  cr c1,c0
  swap c4,c1
  cz c4,c0
  cr c3,c0
  swap c3,c2
  cr c2,c0
  swap c3,c2
  swap c3,c1
  cr c1,c0
  swap c3,c1
  cz c3,c0
  cr c2,c0
  swap c2,c1
  cr c1,c0
  swap c2,c1
  cz c2,c0
  cr c1,c0
  cz c1,c0
  h c0

.qft
  h c0
  cr c1,c0
  cr c2,c0
  cr c3,c0
  cr c4,c0
  cr c5,c0
  cr c6,c0
  h c1
  cr c2,c1
  cr c3,c1
  cr c4,c1
  cr c5,c1
  cr c6,c1
  h c2
  cr c3,c2
  cr c4,c2
  cr c5,c2
  cr c6,c2
  h c3
  cr c4,c3
  cr c5,c3
  cr c6,c3
  h c4
  cr c5,c4
  cr c6,c4
  h c5
  cr c6,c5
  h c6

.cadd
  toffoli a2,b0,z
  cz z,c3
  swap z,c3
  cr c3,c2
  cr c3,c1
  cr c3,c0
  swap z,c3
  toffoli a2,b0,z
  toffoli a2,b1,z
  cz z,c4
  swap z,c4
  cr c4,c3
  cr c4,c2
  cr c4,c1
  cr c4,c0
  swap z,c4
  toffoli a2,b1,z
  toffoli a2,b2,z
  cz z,c5
  swap z,c5
  cr c5,c4
  cr c5,c3
  cr c5,c2
  cr c5,c1
  cr c5,c0
  swap z,c5
  toffoli a2,b2,z
  toffoli a2,b3,z
  cz z,c6
  swap z,c6
  cr c6,c5
  cr c6,c4
  cr c6,c3
  cr c6,c2
  cr c6,c1
  cr c6,c0
  swap z,c6
  toffoli a2,b3,z

.iqft
  h c6
  cr c6,c5
  cz c6,c5
  h c5
  cr c6,c4
  swap c6,c5
  cr c5,c4
  swap c6,c5
  cz c6,c4
  cr c5,c4
  cz c5,c4
  h c4
  cr c6,c3
  swap c6,c5
  cr c5,c3
  swap c6,c5
  swap c6,c4
  cr c4,c3
  swap c6,c4
  cz c6,c3
  cr c5,c3
  swap c5,c4
  cr c4,c3
  swap c5,c4
  cz c5,c3
  cr c4,c3
  cz c4,c3
  h c3
  cr c6,c2
  swap c6,c5
  cr c5,c2
  swap c6,c5
  swap c6,c4
  cr c4,c2
  swap c6,c4
  swap c6,c3
  cr c3,c2
  swap c6,c3
  cz c6,c2
  cr c5,c2
  swap c5,c4
  cr c4,c2
  swap c5,c4
  swap c5,c3
  cr c3,c2
  swap c5,c3
  cz c5,c2
  cr c4,c2
  swap c4,c3
  cr c3,c2
  swap c4,c3
  cz c4,c2
  cr c3,c2
  cz c3,c2
  h c2
  cr c6,c1
  swap c6,c5
  cr c5,c1
  swap c6,c5
  swap c6,c4
  cr c4,c1
  swap c6,c4
  swap c6,c3
  cr c3,c1
  swap c6,c3
  swap c6,c2
  cr c2,c1
  swap c6,c2
  cz c6,c1
  cr c5,c1
  swap c5,c4
  cr c4,c1
  swap c5,c4
  swap c5,c3
  cr c3,c1
  swap c5,c3
  swap c5,c2
  cr c2,c1
  swap c5,c2
  cz c5,c1
  cr c4,c1
  swap c4,c3
  cr c3,c1
  swap c4,c3
  swap c4,c2
  cr c2,c1
  swap c4,c2
  cz c4,c1
  cr c3,c1
  swap c3,c2
  cr c2,c1
  swap c3,c2
  cz c3,c1
  cr c2,c1
  cz c2,c1
  h c1
  cr c6,c0
  swap c6,c5
  cr c5,c0
  swap c6,c5
  swap c6,c4
  cr c4,c0
  swap c6,c4
  swap c6,c3
  cr c3,c0
  swap c6,c3
  swap c6,c2
  cr c2,c0
  swap c6,c2
  swap c6,c1
  cr c1,c0
  swap c6,c1
  cz c6,c0
  cr c5,c0
  swap c5,c4
  cr c4,c0
  swap c5,c4
  swap c5,c3
  cr c3,c0
  swap c5,c3
  swap c5,c2
  cr c2,c0
  swap c5,c2
  swap c5,c1
  cr c1,c0
  swap c5,c1
  cz c5,c0
  cr c4,c0
  swap c4,c3
  cr c3,c0
  swap c4,c3
  swap c4,c2
  cr c2,c0
  swap c4,c2
  swap c4,c1
  cr c1,c0
  swap c4,c1
  cz c4,c0
  cr c3,c0
  swap c3,c2
  cr c2,c0
  swap c3,c2
  swap c3,c1
  cr c1,c0
  swap c3,c1
  cz c3,c0
  cr c2,c0
  swap c2,c1
  cr c1,c0
  swap c2,c1
  cz c2,c0
  cr c1,c0
  cz c1,c0
  h c0

.qft
  h c0
  cr c1,c0
  cr c2,c0
  cr c3,c0
  cr c4,c0
  cr c5,c0
  h c1
  cr c2,c1
  cr c3,c1
  cr c4,c1
  cr c5,c1
  h c2
  cr c3,c2
  cr c4,c2
  cr c5,c2
  h c3
  cr c4,c3
  cr c5,c3
  h c4
  cr c5,c4
  h c5

.cadd
  toffoli a1,b0,z
  cz z,c2
  swap z,c2
  cr c2,c1
  cr c2,c0
  swap z,c2
  toffoli a1,b0,z
  toffoli a1,b1,z
  cz z,c3
  swap z,c3
  cr c3,c2
  cr c3,c1
  cr c3,c0
  swap z,c3
  toffoli a1,b1,z
  toffoli a1,b2,z
  cz z,c4
  swap z,c4
  cr c4,c3
  cr c4,c2
  cr c4,c1
  cr c4,c0
  swap z,c4
  toffoli a1,b2,z
  toffoli a1,b3,z
  cz z,c5
  swap z,c5
  cr c5,c4
  cr c5,c3
  cr c5,c2
  cr c5,c1
  cr c5,c0
  swap z,c5
  toffoli a1,b3,z

.iqft
  h c5
  cr c5,c4
  cz c5,c4
  h c4
  cr c5,c3
  swap c5,c4
  cr c4,c3
  swap c5,c4
  cz c5,c3
  cr c4,c3
  cz c4,c3
  h c3
  cr c5,c2
  swap c5,c4
  cr c4,c2
  swap c5,c4
  swap c5,c3
  cr c3,c2
  swap c5,c3
  cz c5,c2
  cr c4,c2
  swap c4,c3
  cr c3,c2
  swap c4,c3
  cz c4,c2
  cr c3,c2
  cz c3,c2
  h c2
  cr c5,c1
  swap c5,c4
  cr c4,c1
  swap c5,c4
  swap c5,c3
  cr c3,c1
  swap c5,c3
  swap c5,c2
  cr c2,c1
  swap c5,c2
  cz c5,c1
  cr c4,c1
  swap c4,c3
  cr c3,c1
  swap c4,c3
  swap c4,c2
  cr c2,c1
  swap c4,c2
  cz c4,c1
  cr c3,c1
  swap c3,c2
  cr c2,c1
  swap c3,c2
  cz c3,c1
  cr c2,c1
  cz c2,c1
  h c1
  cr c5,c0
  swap c5,c4
  cr c4,c0
  swap c5,c4
  swap c5,c3
  cr c3,c0
  swap c5,c3
  swap c5,c2
  cr c2,c0
  swap c5,c2
  swap c5,c1
  cr c1,c0
  swap c5,c1
  cz c5,c0
  cr c4,c0
  swap c4,c3
  cr c3,c0
  swap c4,c3
  swap c4,c2
  cr c2,c0
  swap c4,c2
  swap c4,c1
  cr c1,c0
  swap c4,c1
  cz c4,c0
  cr c3,c0
  swap c3,c2
  cr c2,c0
  swap c3,c2
  swap c3,c1
  cr c1,c0
  swap c3,c1
  cz c3,c0
  cr c2,c0
  swap c2,c1
  cr c1,c0
  swap c2,c1
  cz c2,c0
  cr c1,c0
  cz c1,c0
  h c0

.qft
  h c0
  cr c1,c0
  cr c2,c0
  cr c3,c0
  cr c4,c0
  h c1
  cr c2,c1
  cr c3,c1
  cr c4,c1
  h c2
  cr c3,c2
  cr c4,c2
  h c3
  cr c4,c3
  h c4

.cadd
  toffoli a0,b0,z
  cz z,c1
  swap z,c1
  cr c1,c0
  swap z,c1
  toffoli a0,b0,z
  toffoli a0,b1,z
  cz z,c2
  swap z,c2
  cr c2,c1
  cr c2,c0
  swap z,c2
  toffoli a0,b1,z
  toffoli a0,b2,z
  cz z,c3
  swap z,c3
  cr c3,c2
  cr c3,c1
  cr c3,c0
  swap z,c3
  toffoli a0,b2,z
  toffoli a0,b3,z
  cz z,c4
  swap z,c4
  cr c4,c3
  cr c4,c2
  cr c4,c1
  cr c4,c0
  swap z,c4
  toffoli a0,b3,z

.iqft
  h c4
  cr c4,c3
  cz c4,c3
  h c3
  cr c4,c2
  swap c4,c3
  cr c3,c2
  swap c4,c3
  cz c4,c2
  cr c3,c2
  cz c3,c2
  h c2
  cr c4,c1
  swap c4,c3
  cr c3,c1
  swap c4,c3
  swap c4,c2
  cr c2,c1
  swap c4,c2
  cz c4,c1
  cr c3,c1
  swap c3,c2
  cr c2,c1
  swap c3,c2
  cz c3,c1
  cr c2,c1
  cz c2,c1
  h c1
  cr c4,c0
  swap c4,c3
  cr c3,c0
  swap c4,c3
  swap c4,c2
  cr c2,c0
  swap c4,c2
  swap c4,c1
  cr c1,c0
  swap c4,c1
  cz c4,c0
  cr c3,c0
  swap c3,c2
  cr c2,c0
  swap c3,c2
  swap c3,c1
  cr c1,c0
  swap c3,c1
  cz c3,c0
  cr c2,c0
  swap c2,c1
  cr c1,c0
  swap c2,c1
  cz c2,c0
  cr c1,c0
  cz c1,c0
  h c0

.result
  measure
  display