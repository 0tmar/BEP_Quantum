#function: Thapliyal Division

qubits 18

.init
  map q0,q0
  map q1,q1
  map q2,q2
  map q3,q3
  map q4,q4
  map q5,q5
  map q6,q6
  map q7,q7
  map q8,q8
  map q9,r0
  map q10,r1
  map q11,r2
  map q12,r3
  map q13,d0
  map q14,d1
  map q15,d2
  map q16,d3
  map q17,c
  x q8
  x d0
  x d1
  x d2
  display

.iteration
  ### Iteration 1, where Y = [q8, ..., r2] and R = r3
  
  # SUB: Y - D, on Y
  x q8
  x r0
  x r1
  x r2
  cx d0,q8
  cx d0,c
  toffoli c,q8,d0
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
  toffoli c,q8,d0
  cx d0,c
  cx c,q8
  x q8
  x r0
  x r1
  x r2
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r3,d0,q8
  cx d0,c
  toffoli c,q8,d0
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
  toffoli c,q8,d0
  cx d0,c
  toffoli r3,c,q8
  
  # x on R
  x r3

.iteration
  ### Iteration 2, where Y = [q7, ..., r1] and R = r2
  
  # SUB: Y - D, on Y
  x q7
  x q8
  x r0
  x r1
  cx d0,q7
  cx d0,c
  toffoli c,q7,d0
  cx d1,q8
  cx d1,d0
  toffoli d0,q8,d1
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
  toffoli d0,q8,d1
  cx d1,d0
  cx d0,q8
  toffoli c,q7,d0
  cx d0,c
  cx c,q7
  x q7
  x q8
  x r0
  x r1
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r2,d0,q7
  cx d0,c
  toffoli c,q7,d0
  toffoli r2,d1,q8
  cx d1,d0
  toffoli d0,q8,d1
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
  toffoli d0,q8,d1
  cx d1,d0
  toffoli r2,d0,q8
  toffoli c,q7,d0
  cx d0,c
  toffoli r2,c,q7
  
  # x on R
  x r2

.iteration
  ### Iteration 3, where Y = [q6, ..., r0] and R = r1
  
  # SUB: Y - D, on Y
  x q6
  x q7
  x q8
  x r0
  cx d0,q6
  cx d0,c
  toffoli c,q6,d0
  cx d1,q7
  cx d1,d0
  toffoli d0,q7,d1
  cx d2,q8
  cx d2,d1
  toffoli d1,q8,d2
  cx d3,r0
  cx d3,d2
  toffoli d2,r0,d3
  cx d3,r1
  toffoli d2,r0,d3
  cx d3,d2
  cx d2,r0
  toffoli d1,q8,d2
  cx d2,d1
  cx d1,q8
  toffoli d0,q7,d1
  cx d1,d0
  cx d0,q7
  toffoli c,q6,d0
  cx d0,c
  cx c,q6
  x q6
  x q7
  x q8
  x r0
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r1,d0,q6
  cx d0,c
  toffoli c,q6,d0
  toffoli r1,d1,q7
  cx d1,d0
  toffoli d0,q7,d1
  toffoli r1,d2,q8
  cx d2,d1
  toffoli d1,q8,d2
  toffoli r1,d3,r0
  cx d3,d2
  toffoli d2,r0,d3
  toffoli d2,r0,d3
  cx d3,d2
  toffoli r1,d2,r0
  toffoli d1,q8,d2
  cx d2,d1
  toffoli r1,d1,q8
  toffoli d0,q7,d1
  cx d1,d0
  toffoli r1,d0,q7
  toffoli c,q6,d0
  cx d0,c
  toffoli r1,c,q6
  
  # x on R
  x r1

.iteration
  ### Iteration 4, where Y = [q5, ..., q8] and R = r0
  
  # SUB: Y - D, on Y
  x q5
  x q6
  x q7
  x q8
  cx d0,q5
  cx d0,c
  toffoli c,q5,d0
  cx d1,q6
  cx d1,d0
  toffoli d0,q6,d1
  cx d2,q7
  cx d2,d1
  toffoli d1,q7,d2
  cx d3,q8
  cx d3,d2
  toffoli d2,q8,d3
  cx d3,r0
  toffoli d2,q8,d3
  cx d3,d2
  cx d2,q8
  toffoli d1,q7,d2
  cx d2,d1
  cx d1,q7
  toffoli d0,q6,d1
  cx d1,d0
  cx d0,q6
  toffoli c,q5,d0
  cx d0,c
  cx c,q5
  x q5
  x q6
  x q7
  x q8
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli r0,d0,q5
  cx d0,c
  toffoli c,q5,d0
  toffoli r0,d1,q6
  cx d1,d0
  toffoli d0,q6,d1
  toffoli r0,d2,q7
  cx d2,d1
  toffoli d1,q7,d2
  toffoli r0,d3,q8
  cx d3,d2
  toffoli d2,q8,d3
  toffoli d2,q8,d3
  cx d3,d2
  toffoli r0,d2,q8
  toffoli d1,q7,d2
  cx d2,d1
  toffoli r0,d1,q7
  toffoli d0,q6,d1
  cx d1,d0
  toffoli r0,d0,q6
  toffoli c,q5,d0
  cx d0,c
  toffoli r0,c,q5
  
  # x on R
  x r0

.iteration
  ### Iteration 5, where Y = [q4, ..., q7] and R = q8
  
  # SUB: Y - D, on Y
  x q4
  x q5
  x q6
  x q7
  cx d0,q4
  cx d0,c
  toffoli c,q4,d0
  cx d1,q5
  cx d1,d0
  toffoli d0,q5,d1
  cx d2,q6
  cx d2,d1
  toffoli d1,q6,d2
  cx d3,q7
  cx d3,d2
  toffoli d2,q7,d3
  cx d3,q8
  toffoli d2,q7,d3
  cx d3,d2
  cx d2,q7
  toffoli d1,q6,d2
  cx d2,d1
  cx d1,q6
  toffoli d0,q5,d1
  cx d1,d0
  cx d0,q5
  toffoli c,q4,d0
  cx d0,c
  cx c,q4
  x q4
  x q5
  x q6
  x q7
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli q8,d0,q4
  cx d0,c
  toffoli c,q4,d0
  toffoli q8,d1,q5
  cx d1,d0
  toffoli d0,q5,d1
  toffoli q8,d2,q6
  cx d2,d1
  toffoli d1,q6,d2
  toffoli q8,d3,q7
  cx d3,d2
  toffoli d2,q7,d3
  toffoli d2,q7,d3
  cx d3,d2
  toffoli q8,d2,q7
  toffoli d1,q6,d2
  cx d2,d1
  toffoli q8,d1,q6
  toffoli d0,q5,d1
  cx d1,d0
  toffoli q8,d0,q5
  toffoli c,q4,d0
  cx d0,c
  toffoli q8,c,q4
  
  # x on R
  x q8

.iteration
  ### Iteration 6, where Y = [q3, ..., q6] and R = q7
  
  # SUB: Y - D, on Y
  x q3
  x q4
  x q5
  x q6
  cx d0,q3
  cx d0,c
  toffoli c,q3,d0
  cx d1,q4
  cx d1,d0
  toffoli d0,q4,d1
  cx d2,q5
  cx d2,d1
  toffoli d1,q5,d2
  cx d3,q6
  cx d3,d2
  toffoli d2,q6,d3
  cx d3,q7
  toffoli d2,q6,d3
  cx d3,d2
  cx d2,q6
  toffoli d1,q5,d2
  cx d2,d1
  cx d1,q5
  toffoli d0,q4,d1
  cx d1,d0
  cx d0,q4
  toffoli c,q3,d0
  cx d0,c
  cx c,q3
  x q3
  x q4
  x q5
  x q6
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli q7,d0,q3
  cx d0,c
  toffoli c,q3,d0
  toffoli q7,d1,q4
  cx d1,d0
  toffoli d0,q4,d1
  toffoli q7,d2,q5
  cx d2,d1
  toffoli d1,q5,d2
  toffoli q7,d3,q6
  cx d3,d2
  toffoli d2,q6,d3
  toffoli d2,q6,d3
  cx d3,d2
  toffoli q7,d2,q6
  toffoli d1,q5,d2
  cx d2,d1
  toffoli q7,d1,q5
  toffoli d0,q4,d1
  cx d1,d0
  toffoli q7,d0,q4
  toffoli c,q3,d0
  cx d0,c
  toffoli q7,c,q3
  
  # x on R
  x q7

.iteration
  ### Iteration 7, where Y = [q2, ..., q5] and R = q6
  
  # SUB: Y - D, on Y
  x q2
  x q3
  x q4
  x q5
  cx d0,q2
  cx d0,c
  toffoli c,q2,d0
  cx d1,q3
  cx d1,d0
  toffoli d0,q3,d1
  cx d2,q4
  cx d2,d1
  toffoli d1,q4,d2
  cx d3,q5
  cx d3,d2
  toffoli d2,q5,d3
  cx d3,q6
  toffoli d2,q5,d3
  cx d3,d2
  cx d2,q5
  toffoli d1,q4,d2
  cx d2,d1
  cx d1,q4
  toffoli d0,q3,d1
  cx d1,d0
  cx d0,q3
  toffoli c,q2,d0
  cx d0,c
  cx c,q2
  x q2
  x q3
  x q4
  x q5
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli q6,d0,q2
  cx d0,c
  toffoli c,q2,d0
  toffoli q6,d1,q3
  cx d1,d0
  toffoli d0,q3,d1
  toffoli q6,d2,q4
  cx d2,d1
  toffoli d1,q4,d2
  toffoli q6,d3,q5
  cx d3,d2
  toffoli d2,q5,d3
  toffoli d2,q5,d3
  cx d3,d2
  toffoli q6,d2,q5
  toffoli d1,q4,d2
  cx d2,d1
  toffoli q6,d1,q4
  toffoli d0,q3,d1
  cx d1,d0
  toffoli q6,d0,q3
  toffoli c,q2,d0
  cx d0,c
  toffoli q6,c,q2
  
  # x on R
  x q6

.iteration
  ### Iteration 8, where Y = [q1, ..., q4] and R = q5
  
  # SUB: Y - D, on Y
  x q1
  x q2
  x q3
  x q4
  cx d0,q1
  cx d0,c
  toffoli c,q1,d0
  cx d1,q2
  cx d1,d0
  toffoli d0,q2,d1
  cx d2,q3
  cx d2,d1
  toffoli d1,q3,d2
  cx d3,q4
  cx d3,d2
  toffoli d2,q4,d3
  cx d3,q5
  toffoli d2,q4,d3
  cx d3,d2
  cx d2,q4
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
  x q4
  
  # cADD: Y (+ D), on Y, ctrl'd by R
  toffoli q5,d0,q1
  cx d0,c
  toffoli c,q1,d0
  toffoli q5,d1,q2
  cx d1,d0
  toffoli d0,q2,d1
  toffoli q5,d2,q3
  cx d2,d1
  toffoli d1,q3,d2
  toffoli q5,d3,q4
  cx d3,d2
  toffoli d2,q4,d3
  toffoli d2,q4,d3
  cx d3,d2
  toffoli q5,d2,q4
  toffoli d1,q3,d2
  cx d2,d1
  toffoli q5,d1,q3
  toffoli d0,q2,d1
  cx d1,d0
  toffoli q5,d0,q2
  toffoli c,q1,d0
  cx d0,c
  toffoli q5,c,q1
  
  # x on R
  x q5

.iteration
  ### Iteration 9, where Y = [q0, ..., q3] and R = q4
  
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
  cx d3,q4
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
  toffoli q4,d0,q0
  cx d0,c
  toffoli c,q0,d0
  toffoli q4,d1,q1
  cx d1,d0
  toffoli d0,q1,d1
  toffoli q4,d2,q2
  cx d2,d1
  toffoli d1,q2,d2
  toffoli q4,d3,q3
  cx d3,d2
  toffoli d2,q3,d3
  toffoli d2,q3,d3
  cx d3,d2
  toffoli q4,d2,q3
  toffoli d1,q2,d2
  cx d2,d1
  toffoli q4,d1,q2
  toffoli d0,q1,d1
  cx d1,d0
  toffoli q4,d0,q1
  toffoli c,q0,d0
  cx d0,c
  toffoli q4,c,q0
  
  # x on R
  x q4

.result
  measure
  display