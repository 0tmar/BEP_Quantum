#function: Thapliyal Division

qubits 13

.init
  map q0,q0
  map q1,q1
  map q2,q2
  map q3,q3
  map q4,r0
  map q5,r1
  map q6,r2
  map q7,r3
  map q8,d0
  map q9,d1
  map q10,d2
  map q11,d3
  map q12,c
  x q0
  x q1
  x q2
  x q3
  x d0
  x d1
  x d2
  x d3
  display

.iteration
  ### Iteration 1, where Y = [q3, ..., r2] and R = r3
  
  # SUB: Y - D, on Y
  x q3
  x r0
  x r1
  x r2
  cx d0,q3
  cx d0,c
  toffoli c,q3,d0
  cx d1,r0
  cx d1,d0
  toffoli d0,r0,d1
  cx d2,r1
  cx d2,d1
  toffoli d1,r1,d2
  cx d3,r2
  cx d3,d2
  toffoli d2,r2,d3
  cx d3,r3
  toffoli d2,r2,d3
  cx d3,d2
  cx d2,r2
  toffoli d1,r1,d2
  cx d2,d1
  cx d1,r1
  toffoli d0,r0,d1
  cx d1,d0
  cx d0,r0
  toffoli c,q3,d0
  cx d0,c
  cx c,q3
  x q3
  x r0
  x r1
  x r2
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r3,d0,q3
  cx d0,c
  toffoli c,q3,d0
  toffoli r3,d1,r0
  cx d1,d0
  toffoli d0,r0,d1
  toffoli r3,d2,r1
  cx d2,d1
  toffoli d1,r1,d2
  toffoli r3,d3,r2
  cx d3,d2
  toffoli d2,r2,d3
  toffoli d2,r2,d3
  cx d3,d2
  toffoli r3,d2,r2
  toffoli d1,r1,d2
  cx d2,d1
  toffoli r3,d1,r1
  toffoli d0,r0,d1
  cx d1,d0
  toffoli r3,d0,r0
  toffoli c,q3,d0
  cx d0,c
  toffoli r3,c,q3
  
  # x on R
  x r3

.iteration
  ### Iteration 2, where Y = [q2, ..., r1] and R = r2
  
  # SUB: Y - D, on Y
  x q2
  x q3
  x r0
  x r1
  cx d0,q2
  cx d0,c
  toffoli c,q2,d0
  cx d1,q3
  cx d1,d0
  toffoli d0,q3,d1
  cx d2,r0
  cx d2,d1
  toffoli d1,r0,d2
  cx d3,r1
  cx d3,d2
  toffoli d2,r1,d3
  cx d3,r2
  toffoli d2,r1,d3
  cx d3,d2
  cx d2,r1
  toffoli d1,r0,d2
  cx d2,d1
  cx d1,r0
  toffoli d0,q3,d1
  cx d1,d0
  cx d0,q3
  toffoli c,q2,d0
  cx d0,c
  cx c,q2
  x q2
  x q3
  x r0
  x r1
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r2,d0,q2
  cx d0,c
  toffoli c,q2,d0
  toffoli r2,d1,q3
  cx d1,d0
  toffoli d0,q3,d1
  toffoli r2,d2,r0
  cx d2,d1
  toffoli d1,r0,d2
  toffoli r2,d3,r1
  cx d3,d2
  toffoli d2,r1,d3
  toffoli d2,r1,d3
  cx d3,d2
  toffoli r2,d2,r1
  toffoli d1,r0,d2
  cx d2,d1
  toffoli r2,d1,r0
  toffoli d0,q3,d1
  cx d1,d0
  toffoli r2,d0,q3
  toffoli c,q2,d0
  cx d0,c
  toffoli r2,c,q2
  
  # x on R
  x r2

.iteration
  ### Iteration 3, where Y = [q1, ..., r0] and R = r1
  
  # SUB: Y - D, on Y
  x q1
  x q2
  x q3
  x r0
  cx d0,q1
  cx d0,c
  toffoli c,q1,d0
  cx d1,q2
  cx d1,d0
  toffoli d0,q2,d1
  cx d2,q3
  cx d2,d1
  toffoli d1,q3,d2
  cx d3,r0
  cx d3,d2
  toffoli d2,r0,d3
  cx d3,r1
  toffoli d2,r0,d3
  cx d3,d2
  cx d2,r0
  toffoli d1,q3,d2
  cx d2,d1
  cx d1,q3
  toffoli d0,q2,d1
  cx d1,d0
  cx d0,q2
  toffoli c,q1,d0
  cx d0,c
  cx c,q1
  x q1
  x q2
  x q3
  x r0
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r1,d0,q1
  cx d0,c
  toffoli c,q1,d0
  toffoli r1,d1,q2
  cx d1,d0
  toffoli d0,q2,d1
  toffoli r1,d2,q3
  cx d2,d1
  toffoli d1,q3,d2
  toffoli r1,d3,r0
  cx d3,d2
  toffoli d2,r0,d3
  toffoli d2,r0,d3
  cx d3,d2
  toffoli r1,d2,r0
  toffoli d1,q3,d2
  cx d2,d1
  toffoli r1,d1,q3
  toffoli d0,q2,d1
  cx d1,d0
  toffoli r1,d0,q2
  toffoli c,q1,d0
  cx d0,c
  toffoli r1,c,q1
  
  # x on R
  x r1

.iteration
  ### Iteration 4, where Y = [q0, ..., q3] and R = r0
  
  # SUB: Y - D, on Y
  x q0
  x q1
  x q2
  x q3
  cx d0,q0
  cx d0,c
  toffoli c,q0,d0
  cx d1,q1
  cx d1,d0
  toffoli d0,q1,d1
  cx d2,q2
  cx d2,d1
  toffoli d1,q2,d2
  cx d3,q3
  cx d3,d2
  toffoli d2,q3,d3
  cx d3,r0
  toffoli d2,q3,d3
  cx d3,d2
  cx d2,q3
  toffoli d1,q2,d2
  cx d2,d1
  cx d1,q2
  toffoli d0,q1,d1
  cx d1,d0
  cx d0,q1
  toffoli c,q0,d0
  cx d0,c
  cx c,q0
  x q0
  x q1
  x q2
  x q3
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r0,d0,q0
  cx d0,c
  toffoli c,q0,d0
  toffoli r0,d1,q1
  cx d1,d0
  toffoli d0,q1,d1
  toffoli r0,d2,q2
  cx d2,d1
  toffoli d1,q2,d2
  toffoli r0,d3,q3
  cx d3,d2
  toffoli d2,q3,d3
  toffoli d2,q3,d3
  cx d3,d2
  toffoli r0,d2,q3
  toffoli d1,q2,d2
  cx d2,d1
  toffoli r0,d1,q2
  toffoli d0,q1,d1
  cx d1,d0
  toffoli r0,d0,q1
  toffoli c,q0,d0
  cx d0,c
  toffoli r0,c,q0
  
  # x on R
  x r0

.result
  measure
  display