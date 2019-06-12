#function: Thapliyal Division

qubits 14

.init
  map q0,n0
  map q1,n1
  map q2,n2
  map q3,n3
  map q4,n4
  map q5,n5
  map q6,o0
  map q7,o1
  map q8,o2
  map q9,o3
  map q10,d0
  map q11,d1
  map q12,d2
  map q13,d3
  x n5
  x d1
  display

.thapliyal_division
  ##### Thapliyal division with N = [n0, ..., n5], O = [o0, ..., o3] and D = [d0, ..., d3]
  
  ### Iteration 1, where Y = [n5, ..., o2] and Z = o3
  
  # SUB: Y - D, on Y
  x n5
  x o0
  x o1
  x o2
  cx d1,o0
  cx d2,o1
  cx d3,o2
  cx d3,o3
  cx d2,d3
  cx d1,d2
  toffoli d0,n5,d1
  toffoli d1,o0,d2
  toffoli d2,o1,d3
  toffoli d3,o2,o3
  cx d3,o2
  toffoli d2,o1,d3
  cx d2,o1
  toffoli d1,o0,d2
  cx d1,o0
  toffoli d0,n5,d1
  cx d0,n5
  cx d1,d2
  cx d2,d3
  cx d3,o2
  cx d2,o1
  cx d1,o0
  x n5
  x o0
  x o1
  x o2
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,o0
  cx d2,o1
  cx d3,o2
  cx d2,d3
  cx d1,d2
  toffoli d0,n5,d1
  toffoli d1,o0,d2
  toffoli d2,o1,d3
  toffoli o3,d3,o2
  toffoli d2,o1,d3
  toffoli o3,d2,o1
  toffoli d1,o0,d2
  toffoli o3,d1,o0
  toffoli d0,n5,d1
  toffoli o3,d0,n5
  cx d1,d2
  cx d2,d3
  cx d3,o2
  cx d2,o1
  cx d1,o0
  
  # x on Z
  x o3
  
  ### Iteration 2, where Y = [n4, ..., o1] and Z = o2
  
  # SUB: Y - D, on Y
  x n4
  x n5
  x o0
  x o1
  cx d1,n5
  cx d2,o0
  cx d3,o1
  cx d3,o2
  cx d2,d3
  cx d1,d2
  toffoli d0,n4,d1
  toffoli d1,n5,d2
  toffoli d2,o0,d3
  toffoli d3,o1,o2
  cx d3,o1
  toffoli d2,o0,d3
  cx d2,o0
  toffoli d1,n5,d2
  cx d1,n5
  toffoli d0,n4,d1
  cx d0,n4
  cx d1,d2
  cx d2,d3
  cx d3,o1
  cx d2,o0
  cx d1,n5
  x n4
  x n5
  x o0
  x o1
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n5
  cx d2,o0
  cx d3,o1
  cx d2,d3
  cx d1,d2
  toffoli d0,n4,d1
  toffoli d1,n5,d2
  toffoli d2,o0,d3
  toffoli o2,d3,o1
  toffoli d2,o0,d3
  toffoli o2,d2,o0
  toffoli d1,n5,d2
  toffoli o2,d1,n5
  toffoli d0,n4,d1
  toffoli o2,d0,n4
  cx d1,d2
  cx d2,d3
  cx d3,o1
  cx d2,o0
  cx d1,n5
  
  # x on Z
  x o2
  
  ### Iteration 3, where Y = [n3, ..., o0] and Z = o1
  
  # SUB: Y - D, on Y
  x n3
  x n4
  x n5
  x o0
  cx d1,n4
  cx d2,n5
  cx d3,o0
  cx d3,o1
  cx d2,d3
  cx d1,d2
  toffoli d0,n3,d1
  toffoli d1,n4,d2
  toffoli d2,n5,d3
  toffoli d3,o0,o1
  cx d3,o0
  toffoli d2,n5,d3
  cx d2,n5
  toffoli d1,n4,d2
  cx d1,n4
  toffoli d0,n3,d1
  cx d0,n3
  cx d1,d2
  cx d2,d3
  cx d3,o0
  cx d2,n5
  cx d1,n4
  x n3
  x n4
  x n5
  x o0
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n4
  cx d2,n5
  cx d3,o0
  cx d2,d3
  cx d1,d2
  toffoli d0,n3,d1
  toffoli d1,n4,d2
  toffoli d2,n5,d3
  toffoli o1,d3,o0
  toffoli d2,n5,d3
  toffoli o1,d2,n5
  toffoli d1,n4,d2
  toffoli o1,d1,n4
  toffoli d0,n3,d1
  toffoli o1,d0,n3
  cx d1,d2
  cx d2,d3
  cx d3,o0
  cx d2,n5
  cx d1,n4
  
  # x on Z
  x o1
  
  ### Iteration 4, where Y = [n2, ..., n5] and Z = o0
  
  # SUB: Y - D, on Y
  x n2
  x n3
  x n4
  x n5
  cx d1,n3
  cx d2,n4
  cx d3,n5
  cx d3,o0
  cx d2,d3
  cx d1,d2
  toffoli d0,n2,d1
  toffoli d1,n3,d2
  toffoli d2,n4,d3
  toffoli d3,n5,o0
  cx d3,n5
  toffoli d2,n4,d3
  cx d2,n4
  toffoli d1,n3,d2
  cx d1,n3
  toffoli d0,n2,d1
  cx d0,n2
  cx d1,d2
  cx d2,d3
  cx d3,n5
  cx d2,n4
  cx d1,n3
  x n2
  x n3
  x n4
  x n5
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n3
  cx d2,n4
  cx d3,n5
  cx d2,d3
  cx d1,d2
  toffoli d0,n2,d1
  toffoli d1,n3,d2
  toffoli d2,n4,d3
  toffoli o0,d3,n5
  toffoli d2,n4,d3
  toffoli o0,d2,n4
  toffoli d1,n3,d2
  toffoli o0,d1,n3
  toffoli d0,n2,d1
  toffoli o0,d0,n2
  cx d1,d2
  cx d2,d3
  cx d3,n5
  cx d2,n4
  cx d1,n3
  
  # x on Z
  x o0
  
  ### Iteration 5, where Y = [n1, ..., n4] and Z = n5
  
  # SUB: Y - D, on Y
  x n1
  x n2
  x n3
  x n4
  cx d1,n2
  cx d2,n3
  cx d3,n4
  cx d3,n5
  cx d2,d3
  cx d1,d2
  toffoli d0,n1,d1
  toffoli d1,n2,d2
  toffoli d2,n3,d3
  toffoli d3,n4,n5
  cx d3,n4
  toffoli d2,n3,d3
  cx d2,n3
  toffoli d1,n2,d2
  cx d1,n2
  toffoli d0,n1,d1
  cx d0,n1
  cx d1,d2
  cx d2,d3
  cx d3,n4
  cx d2,n3
  cx d1,n2
  x n1
  x n2
  x n3
  x n4
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n2
  cx d2,n3
  cx d3,n4
  cx d2,d3
  cx d1,d2
  toffoli d0,n1,d1
  toffoli d1,n2,d2
  toffoli d2,n3,d3
  toffoli n5,d3,n4
  toffoli d2,n3,d3
  toffoli n5,d2,n3
  toffoli d1,n2,d2
  toffoli n5,d1,n2
  toffoli d0,n1,d1
  toffoli n5,d0,n1
  cx d1,d2
  cx d2,d3
  cx d3,n4
  cx d2,n3
  cx d1,n2
  
  # x on Z
  x n5
  
  ### Iteration 6, where Y = [n0, ..., n3] and Z = n4
  
  # SUB: Y - D, on Y
  x n0
  x n1
  x n2
  x n3
  cx d1,n1
  cx d2,n2
  cx d3,n3
  cx d3,n4
  cx d2,d3
  cx d1,d2
  toffoli d0,n0,d1
  toffoli d1,n1,d2
  toffoli d2,n2,d3
  toffoli d3,n3,n4
  cx d3,n3
  toffoli d2,n2,d3
  cx d2,n2
  toffoli d1,n1,d2
  cx d1,n1
  toffoli d0,n0,d1
  cx d0,n0
  cx d1,d2
  cx d2,d3
  cx d3,n3
  cx d2,n2
  cx d1,n1
  x n0
  x n1
  x n2
  x n3
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n1
  cx d2,n2
  cx d3,n3
  cx d2,d3
  cx d1,d2
  toffoli d0,n0,d1
  toffoli d1,n1,d2
  toffoli d2,n2,d3
  toffoli n4,d3,n3
  toffoli d2,n2,d3
  toffoli n4,d2,n2
  toffoli d1,n1,d2
  toffoli n4,d1,n1
  toffoli d0,n0,d1
  toffoli n4,d0,n0
  cx d1,d2
  cx d2,d3
  cx d3,n3
  cx d2,n2
  cx d1,n1
  
  # x on Z
  x n4

.result
  measure
  display