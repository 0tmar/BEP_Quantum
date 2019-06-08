#function: Thapliyal Division

qubits 14

.init
  map q0,q0
  map q1,q1
  map q2,q2
  map q3,q3
  map q4,q4
  map q5,q5
  map q6,r0
  map q7,r1
  map q8,r2
  map q9,r3
  map q10,d0
  map q11,d1
  map q12,d2
  map q13,d3
  x q5
  x d1
  display

.thapliyal_division
  ##### Thapliyal division with Q = [q0, ..., q5], R = [r0, ..., r3] and D = [d0, ..., d3]
  
  ### Iteration 1, where Y = [q5, ..., r2] and Z = r3
  
  # SUB: Y - D, on Y
  x q5
  x r0
  x r1
  x r2
  cx d1,r0
  cx d2,r1
  cx d3,r2
  cx d3,r3
  cx d2,d3
  cx d1,d2
  toffoli d0,q5,d1
  toffoli d1,r0,d2
  toffoli d2,r1,d3
  toffoli d3,r2,r3
  cx d3,r2
  toffoli d2,r1,d3
  cx d2,r1
  toffoli d1,r0,d2
  cx d1,r0
  toffoli d0,q5,d1
  cx d0,q5
  cx d1,d2
  cx d2,d3
  cx d3,r2
  cx d2,r1
  cx d1,r0
  x q5
  x r0
  x r1
  x r2
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,r0
  cx d2,r1
  cx d3,r2
  cx d2,d3
  cx d1,d2
  toffoli d0,q5,d1
  toffoli d1,r0,d2
  toffoli d2,r1,d3
  toffoli r3,d3,r2
  toffoli d2,r1,d3
  toffoli r3,d2,r1
  toffoli d1,r0,d2
  toffoli r3,d1,r0
  toffoli d0,q5,d1
  toffoli r3,d0,q5
  cx d1,d2
  cx d2,d3
  cx d3,r2
  cx d2,r1
  cx d1,r0
  
  # x on Z
  x r3
  
  ### Iteration 2, where Y = [q4, ..., r1] and Z = r2
  
  # SUB: Y - D, on Y
  x q4
  x q5
  x r0
  x r1
  cx d1,q5
  cx d2,r0
  cx d3,r1
  cx d3,r2
  cx d2,d3
  cx d1,d2
  toffoli d0,q4,d1
  toffoli d1,q5,d2
  toffoli d2,r0,d3
  toffoli d3,r1,r2
  cx d3,r1
  toffoli d2,r0,d3
  cx d2,r0
  toffoli d1,q5,d2
  cx d1,q5
  toffoli d0,q4,d1
  cx d0,q4
  cx d1,d2
  cx d2,d3
  cx d3,r1
  cx d2,r0
  cx d1,q5
  x q4
  x q5
  x r0
  x r1
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,q5
  cx d2,r0
  cx d3,r1
  cx d2,d3
  cx d1,d2
  toffoli d0,q4,d1
  toffoli d1,q5,d2
  toffoli d2,r0,d3
  toffoli r2,d3,r1
  toffoli d2,r0,d3
  toffoli r2,d2,r0
  toffoli d1,q5,d2
  toffoli r2,d1,q5
  toffoli d0,q4,d1
  toffoli r2,d0,q4
  cx d1,d2
  cx d2,d3
  cx d3,r1
  cx d2,r0
  cx d1,q5
  
  # x on Z
  x r2
  
  ### Iteration 3, where Y = [q3, ..., r0] and Z = r1
  
  # SUB: Y - D, on Y
  x q3
  x q4
  x q5
  x r0
  cx d1,q4
  cx d2,q5
  cx d3,r0
  cx d3,r1
  cx d2,d3
  cx d1,d2
  toffoli d0,q3,d1
  toffoli d1,q4,d2
  toffoli d2,q5,d3
  toffoli d3,r0,r1
  cx d3,r0
  toffoli d2,q5,d3
  cx d2,q5
  toffoli d1,q4,d2
  cx d1,q4
  toffoli d0,q3,d1
  cx d0,q3
  cx d1,d2
  cx d2,d3
  cx d3,r0
  cx d2,q5
  cx d1,q4
  x q3
  x q4
  x q5
  x r0
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,q4
  cx d2,q5
  cx d3,r0
  cx d2,d3
  cx d1,d2
  toffoli d0,q3,d1
  toffoli d1,q4,d2
  toffoli d2,q5,d3
  toffoli r1,d3,r0
  toffoli d2,q5,d3
  toffoli r1,d2,q5
  toffoli d1,q4,d2
  toffoli r1,d1,q4
  toffoli d0,q3,d1
  toffoli r1,d0,q3
  cx d1,d2
  cx d2,d3
  cx d3,r0
  cx d2,q5
  cx d1,q4
  
  # x on Z
  x r1
  
  ### Iteration 4, where Y = [q2, ..., q5] and Z = r0
  
  # SUB: Y - D, on Y
  x q2
  x q3
  x q4
  x q5
  cx d1,q3
  cx d2,q4
  cx d3,q5
  cx d3,r0
  cx d2,d3
  cx d1,d2
  toffoli d0,q2,d1
  toffoli d1,q3,d2
  toffoli d2,q4,d3
  toffoli d3,q5,r0
  cx d3,q5
  toffoli d2,q4,d3
  cx d2,q4
  toffoli d1,q3,d2
  cx d1,q3
  toffoli d0,q2,d1
  cx d0,q2
  cx d1,d2
  cx d2,d3
  cx d3,q5
  cx d2,q4
  cx d1,q3
  x q2
  x q3
  x q4
  x q5
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,q3
  cx d2,q4
  cx d3,q5
  cx d2,d3
  cx d1,d2
  toffoli d0,q2,d1
  toffoli d1,q3,d2
  toffoli d2,q4,d3
  toffoli r0,d3,q5
  toffoli d2,q4,d3
  toffoli r0,d2,q4
  toffoli d1,q3,d2
  toffoli r0,d1,q3
  toffoli d0,q2,d1
  toffoli r0,d0,q2
  cx d1,d2
  cx d2,d3
  cx d3,q5
  cx d2,q4
  cx d1,q3
  
  # x on Z
  x r0
  
  ### Iteration 5, where Y = [q1, ..., q4] and Z = q5
  
  # SUB: Y - D, on Y
  x q1
  x q2
  x q3
  x q4
  cx d1,q2
  cx d2,q3
  cx d3,q4
  cx d3,q5
  cx d2,d3
  cx d1,d2
  toffoli d0,q1,d1
  toffoli d1,q2,d2
  toffoli d2,q3,d3
  toffoli d3,q4,q5
  cx d3,q4
  toffoli d2,q3,d3
  cx d2,q3
  toffoli d1,q2,d2
  cx d1,q2
  toffoli d0,q1,d1
  cx d0,q1
  cx d1,d2
  cx d2,d3
  cx d3,q4
  cx d2,q3
  cx d1,q2
  x q1
  x q2
  x q3
  x q4
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,q2
  cx d2,q3
  cx d3,q4
  cx d2,d3
  cx d1,d2
  toffoli d0,q1,d1
  toffoli d1,q2,d2
  toffoli d2,q3,d3
  toffoli q5,d3,q4
  toffoli d2,q3,d3
  toffoli q5,d2,q3
  toffoli d1,q2,d2
  toffoli q5,d1,q2
  toffoli d0,q1,d1
  toffoli q5,d0,q1
  cx d1,d2
  cx d2,d3
  cx d3,q4
  cx d2,q3
  cx d1,q2
  
  # x on Z
  x q5
  
  ### Iteration 6, where Y = [q0, ..., q3] and Z = q4
  
  # SUB: Y - D, on Y
  x q0
  x q1
  x q2
  x q3
  cx d1,q1
  cx d2,q2
  cx d3,q3
  cx d3,q4
  cx d2,d3
  cx d1,d2
  toffoli d0,q0,d1
  toffoli d1,q1,d2
  toffoli d2,q2,d3
  toffoli d3,q3,q4
  cx d3,q3
  toffoli d2,q2,d3
  cx d2,q2
  toffoli d1,q1,d2
  cx d1,q1
  toffoli d0,q0,d1
  cx d0,q0
  cx d1,d2
  cx d2,d3
  cx d3,q3
  cx d2,q2
  cx d1,q1
  x q0
  x q1
  x q2
  x q3
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,q1
  cx d2,q2
  cx d3,q3
  cx d2,d3
  cx d1,d2
  toffoli d0,q0,d1
  toffoli d1,q1,d2
  toffoli d2,q2,d3
  toffoli q4,d3,q3
  toffoli d2,q2,d3
  toffoli q4,d2,q2
  toffoli d1,q1,d2
  toffoli q4,d1,q1
  toffoli d0,q0,d1
  toffoli q4,d0,q0
  cx d1,d2
  cx d2,d3
  cx d3,q3
  cx d2,q2
  cx d1,q1
  
  # x on Z
  x q4

.result
  measure
  display