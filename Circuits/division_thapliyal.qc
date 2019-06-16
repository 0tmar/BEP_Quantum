#function: Thapliyal Division

qubits 25

.init
  map q0,n0
  map q1,n1
  map q2,n2
  map q3,n3
  map q4,n4
  map q5,n5
  map q6,n6
  map q7,n7
  map q8,n8
  map q9,n9
  map q10,n10
  map q11,n11
  map q12,n12
  map q13,n13
  map q14,n14
  map q15,o0
  map q16,o1
  map q17,o2
  map q18,o3
  map q19,o4
  map q20,d0
  map q21,d1
  map q22,d2
  map q23,d3
  map q24,d4
  x n14
  x d0
  x d1
  x d3
  display

.thapliyal_division
  ##### Thapliyal division with N = [n0, ..., n14], O = [o0, ..., o4] and D = [d0, ..., d4]
  
  ### Iteration 1, where Y = [n14, ..., o3] and Z = o4
  
  # SUB: Y - D, on Y
  x n14
  x o0
  x o1
  x o2
  x o3
  cx d1,o0
  cx d2,o1
  cx d3,o2
  cx d4,o3
  cx d4,o4
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n14,d1
  toffoli d1,o0,d2
  toffoli d2,o1,d3
  toffoli d3,o2,d4
  toffoli d4,o3,o4
  cx d4,o3
  toffoli d3,o2,d4
  cx d3,o2
  toffoli d2,o1,d3
  cx d2,o1
  toffoli d1,o0,d2
  cx d1,o0
  toffoli d0,n14,d1
  cx d0,n14
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o3
  cx d3,o2
  cx d2,o1
  cx d1,o0
  x n14
  x o0
  x o1
  x o2
  x o3
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,o0
  cx d2,o1
  cx d3,o2
  cx d4,o3
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n14,d1
  toffoli d1,o0,d2
  toffoli d2,o1,d3
  toffoli d3,o2,d4
  toffoli o4,d4,o3
  toffoli d3,o2,d4
  toffoli o4,d3,o2
  toffoli d2,o1,d3
  toffoli o4,d2,o1
  toffoli d1,o0,d2
  toffoli o4,d1,o0
  toffoli d0,n14,d1
  toffoli o4,d0,n14
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o3
  cx d3,o2
  cx d2,o1
  cx d1,o0
  
  # x on Z
  x o4
  
  ### Iteration 2, where Y = [n13, ..., o2] and Z = o3
  
  # SUB: Y - D, on Y
  x n13
  x n14
  x o0
  x o1
  x o2
  cx d1,n14
  cx d2,o0
  cx d3,o1
  cx d4,o2
  cx d4,o3
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n13,d1
  toffoli d1,n14,d2
  toffoli d2,o0,d3
  toffoli d3,o1,d4
  toffoli d4,o2,o3
  cx d4,o2
  toffoli d3,o1,d4
  cx d3,o1
  toffoli d2,o0,d3
  cx d2,o0
  toffoli d1,n14,d2
  cx d1,n14
  toffoli d0,n13,d1
  cx d0,n13
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o2
  cx d3,o1
  cx d2,o0
  cx d1,n14
  x n13
  x n14
  x o0
  x o1
  x o2
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n14
  cx d2,o0
  cx d3,o1
  cx d4,o2
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n13,d1
  toffoli d1,n14,d2
  toffoli d2,o0,d3
  toffoli d3,o1,d4
  toffoli o3,d4,o2
  toffoli d3,o1,d4
  toffoli o3,d3,o1
  toffoli d2,o0,d3
  toffoli o3,d2,o0
  toffoli d1,n14,d2
  toffoli o3,d1,n14
  toffoli d0,n13,d1
  toffoli o3,d0,n13
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o2
  cx d3,o1
  cx d2,o0
  cx d1,n14
  
  # x on Z
  x o3
  
  ### Iteration 3, where Y = [n12, ..., o1] and Z = o2
  
  # SUB: Y - D, on Y
  x n12
  x n13
  x n14
  x o0
  x o1
  cx d1,n13
  cx d2,n14
  cx d3,o0
  cx d4,o1
  cx d4,o2
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n12,d1
  toffoli d1,n13,d2
  toffoli d2,n14,d3
  toffoli d3,o0,d4
  toffoli d4,o1,o2
  cx d4,o1
  toffoli d3,o0,d4
  cx d3,o0
  toffoli d2,n14,d3
  cx d2,n14
  toffoli d1,n13,d2
  cx d1,n13
  toffoli d0,n12,d1
  cx d0,n12
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o1
  cx d3,o0
  cx d2,n14
  cx d1,n13
  x n12
  x n13
  x n14
  x o0
  x o1
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n13
  cx d2,n14
  cx d3,o0
  cx d4,o1
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n12,d1
  toffoli d1,n13,d2
  toffoli d2,n14,d3
  toffoli d3,o0,d4
  toffoli o2,d4,o1
  toffoli d3,o0,d4
  toffoli o2,d3,o0
  toffoli d2,n14,d3
  toffoli o2,d2,n14
  toffoli d1,n13,d2
  toffoli o2,d1,n13
  toffoli d0,n12,d1
  toffoli o2,d0,n12
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o1
  cx d3,o0
  cx d2,n14
  cx d1,n13
  
  # x on Z
  x o2
  
  ### Iteration 4, where Y = [n11, ..., o0] and Z = o1
  
  # SUB: Y - D, on Y
  x n11
  x n12
  x n13
  x n14
  x o0
  cx d1,n12
  cx d2,n13
  cx d3,n14
  cx d4,o0
  cx d4,o1
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n11,d1
  toffoli d1,n12,d2
  toffoli d2,n13,d3
  toffoli d3,n14,d4
  toffoli d4,o0,o1
  cx d4,o0
  toffoli d3,n14,d4
  cx d3,n14
  toffoli d2,n13,d3
  cx d2,n13
  toffoli d1,n12,d2
  cx d1,n12
  toffoli d0,n11,d1
  cx d0,n11
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o0
  cx d3,n14
  cx d2,n13
  cx d1,n12
  x n11
  x n12
  x n13
  x n14
  x o0
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n12
  cx d2,n13
  cx d3,n14
  cx d4,o0
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n11,d1
  toffoli d1,n12,d2
  toffoli d2,n13,d3
  toffoli d3,n14,d4
  toffoli o1,d4,o0
  toffoli d3,n14,d4
  toffoli o1,d3,n14
  toffoli d2,n13,d3
  toffoli o1,d2,n13
  toffoli d1,n12,d2
  toffoli o1,d1,n12
  toffoli d0,n11,d1
  toffoli o1,d0,n11
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,o0
  cx d3,n14
  cx d2,n13
  cx d1,n12
  
  # x on Z
  x o1
  
  ### Iteration 5, where Y = [n10, ..., n14] and Z = o0
  
  # SUB: Y - D, on Y
  x n10
  x n11
  x n12
  x n13
  x n14
  cx d1,n11
  cx d2,n12
  cx d3,n13
  cx d4,n14
  cx d4,o0
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n10,d1
  toffoli d1,n11,d2
  toffoli d2,n12,d3
  toffoli d3,n13,d4
  toffoli d4,n14,o0
  cx d4,n14
  toffoli d3,n13,d4
  cx d3,n13
  toffoli d2,n12,d3
  cx d2,n12
  toffoli d1,n11,d2
  cx d1,n11
  toffoli d0,n10,d1
  cx d0,n10
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n14
  cx d3,n13
  cx d2,n12
  cx d1,n11
  x n10
  x n11
  x n12
  x n13
  x n14
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n11
  cx d2,n12
  cx d3,n13
  cx d4,n14
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n10,d1
  toffoli d1,n11,d2
  toffoli d2,n12,d3
  toffoli d3,n13,d4
  toffoli o0,d4,n14
  toffoli d3,n13,d4
  toffoli o0,d3,n13
  toffoli d2,n12,d3
  toffoli o0,d2,n12
  toffoli d1,n11,d2
  toffoli o0,d1,n11
  toffoli d0,n10,d1
  toffoli o0,d0,n10
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n14
  cx d3,n13
  cx d2,n12
  cx d1,n11
  
  # x on Z
  x o0
  
  ### Iteration 6, where Y = [n9, ..., n13] and Z = n14
  
  # SUB: Y - D, on Y
  x n9
  x n10
  x n11
  x n12
  x n13
  cx d1,n10
  cx d2,n11
  cx d3,n12
  cx d4,n13
  cx d4,n14
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n9,d1
  toffoli d1,n10,d2
  toffoli d2,n11,d3
  toffoli d3,n12,d4
  toffoli d4,n13,n14
  cx d4,n13
  toffoli d3,n12,d4
  cx d3,n12
  toffoli d2,n11,d3
  cx d2,n11
  toffoli d1,n10,d2
  cx d1,n10
  toffoli d0,n9,d1
  cx d0,n9
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n13
  cx d3,n12
  cx d2,n11
  cx d1,n10
  x n9
  x n10
  x n11
  x n12
  x n13
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n10
  cx d2,n11
  cx d3,n12
  cx d4,n13
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n9,d1
  toffoli d1,n10,d2
  toffoli d2,n11,d3
  toffoli d3,n12,d4
  toffoli n14,d4,n13
  toffoli d3,n12,d4
  toffoli n14,d3,n12
  toffoli d2,n11,d3
  toffoli n14,d2,n11
  toffoli d1,n10,d2
  toffoli n14,d1,n10
  toffoli d0,n9,d1
  toffoli n14,d0,n9
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n13
  cx d3,n12
  cx d2,n11
  cx d1,n10
  
  # x on Z
  x n14
  
  ### Iteration 7, where Y = [n8, ..., n12] and Z = n13
  
  # SUB: Y - D, on Y
  x n8
  x n9
  x n10
  x n11
  x n12
  cx d1,n9
  cx d2,n10
  cx d3,n11
  cx d4,n12
  cx d4,n13
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n8,d1
  toffoli d1,n9,d2
  toffoli d2,n10,d3
  toffoli d3,n11,d4
  toffoli d4,n12,n13
  cx d4,n12
  toffoli d3,n11,d4
  cx d3,n11
  toffoli d2,n10,d3
  cx d2,n10
  toffoli d1,n9,d2
  cx d1,n9
  toffoli d0,n8,d1
  cx d0,n8
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n12
  cx d3,n11
  cx d2,n10
  cx d1,n9
  x n8
  x n9
  x n10
  x n11
  x n12
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n9
  cx d2,n10
  cx d3,n11
  cx d4,n12
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n8,d1
  toffoli d1,n9,d2
  toffoli d2,n10,d3
  toffoli d3,n11,d4
  toffoli n13,d4,n12
  toffoli d3,n11,d4
  toffoli n13,d3,n11
  toffoli d2,n10,d3
  toffoli n13,d2,n10
  toffoli d1,n9,d2
  toffoli n13,d1,n9
  toffoli d0,n8,d1
  toffoli n13,d0,n8
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n12
  cx d3,n11
  cx d2,n10
  cx d1,n9
  
  # x on Z
  x n13
  
  ### Iteration 8, where Y = [n7, ..., n11] and Z = n12
  
  # SUB: Y - D, on Y
  x n7
  x n8
  x n9
  x n10
  x n11
  cx d1,n8
  cx d2,n9
  cx d3,n10
  cx d4,n11
  cx d4,n12
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n7,d1
  toffoli d1,n8,d2
  toffoli d2,n9,d3
  toffoli d3,n10,d4
  toffoli d4,n11,n12
  cx d4,n11
  toffoli d3,n10,d4
  cx d3,n10
  toffoli d2,n9,d3
  cx d2,n9
  toffoli d1,n8,d2
  cx d1,n8
  toffoli d0,n7,d1
  cx d0,n7
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n11
  cx d3,n10
  cx d2,n9
  cx d1,n8
  x n7
  x n8
  x n9
  x n10
  x n11
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n8
  cx d2,n9
  cx d3,n10
  cx d4,n11
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n7,d1
  toffoli d1,n8,d2
  toffoli d2,n9,d3
  toffoli d3,n10,d4
  toffoli n12,d4,n11
  toffoli d3,n10,d4
  toffoli n12,d3,n10
  toffoli d2,n9,d3
  toffoli n12,d2,n9
  toffoli d1,n8,d2
  toffoli n12,d1,n8
  toffoli d0,n7,d1
  toffoli n12,d0,n7
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n11
  cx d3,n10
  cx d2,n9
  cx d1,n8
  
  # x on Z
  x n12
  
  ### Iteration 9, where Y = [n6, ..., n10] and Z = n11
  
  # SUB: Y - D, on Y
  x n6
  x n7
  x n8
  x n9
  x n10
  cx d1,n7
  cx d2,n8
  cx d3,n9
  cx d4,n10
  cx d4,n11
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n6,d1
  toffoli d1,n7,d2
  toffoli d2,n8,d3
  toffoli d3,n9,d4
  toffoli d4,n10,n11
  cx d4,n10
  toffoli d3,n9,d4
  cx d3,n9
  toffoli d2,n8,d3
  cx d2,n8
  toffoli d1,n7,d2
  cx d1,n7
  toffoli d0,n6,d1
  cx d0,n6
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n10
  cx d3,n9
  cx d2,n8
  cx d1,n7
  x n6
  x n7
  x n8
  x n9
  x n10
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n7
  cx d2,n8
  cx d3,n9
  cx d4,n10
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n6,d1
  toffoli d1,n7,d2
  toffoli d2,n8,d3
  toffoli d3,n9,d4
  toffoli n11,d4,n10
  toffoli d3,n9,d4
  toffoli n11,d3,n9
  toffoli d2,n8,d3
  toffoli n11,d2,n8
  toffoli d1,n7,d2
  toffoli n11,d1,n7
  toffoli d0,n6,d1
  toffoli n11,d0,n6
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n10
  cx d3,n9
  cx d2,n8
  cx d1,n7
  
  # x on Z
  x n11
  
  ### Iteration 10, where Y = [n5, ..., n9] and Z = n10
  
  # SUB: Y - D, on Y
  x n5
  x n6
  x n7
  x n8
  x n9
  cx d1,n6
  cx d2,n7
  cx d3,n8
  cx d4,n9
  cx d4,n10
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n5,d1
  toffoli d1,n6,d2
  toffoli d2,n7,d3
  toffoli d3,n8,d4
  toffoli d4,n9,n10
  cx d4,n9
  toffoli d3,n8,d4
  cx d3,n8
  toffoli d2,n7,d3
  cx d2,n7
  toffoli d1,n6,d2
  cx d1,n6
  toffoli d0,n5,d1
  cx d0,n5
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n9
  cx d3,n8
  cx d2,n7
  cx d1,n6
  x n5
  x n6
  x n7
  x n8
  x n9
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n6
  cx d2,n7
  cx d3,n8
  cx d4,n9
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n5,d1
  toffoli d1,n6,d2
  toffoli d2,n7,d3
  toffoli d3,n8,d4
  toffoli n10,d4,n9
  toffoli d3,n8,d4
  toffoli n10,d3,n8
  toffoli d2,n7,d3
  toffoli n10,d2,n7
  toffoli d1,n6,d2
  toffoli n10,d1,n6
  toffoli d0,n5,d1
  toffoli n10,d0,n5
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n9
  cx d3,n8
  cx d2,n7
  cx d1,n6
  
  # x on Z
  x n10
  
  ### Iteration 11, where Y = [n4, ..., n8] and Z = n9
  
  # SUB: Y - D, on Y
  x n4
  x n5
  x n6
  x n7
  x n8
  cx d1,n5
  cx d2,n6
  cx d3,n7
  cx d4,n8
  cx d4,n9
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n4,d1
  toffoli d1,n5,d2
  toffoli d2,n6,d3
  toffoli d3,n7,d4
  toffoli d4,n8,n9
  cx d4,n8
  toffoli d3,n7,d4
  cx d3,n7
  toffoli d2,n6,d3
  cx d2,n6
  toffoli d1,n5,d2
  cx d1,n5
  toffoli d0,n4,d1
  cx d0,n4
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n8
  cx d3,n7
  cx d2,n6
  cx d1,n5
  x n4
  x n5
  x n6
  x n7
  x n8
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n5
  cx d2,n6
  cx d3,n7
  cx d4,n8
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n4,d1
  toffoli d1,n5,d2
  toffoli d2,n6,d3
  toffoli d3,n7,d4
  toffoli n9,d4,n8
  toffoli d3,n7,d4
  toffoli n9,d3,n7
  toffoli d2,n6,d3
  toffoli n9,d2,n6
  toffoli d1,n5,d2
  toffoli n9,d1,n5
  toffoli d0,n4,d1
  toffoli n9,d0,n4
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n8
  cx d3,n7
  cx d2,n6
  cx d1,n5
  
  # x on Z
  x n9
  
  ### Iteration 12, where Y = [n3, ..., n7] and Z = n8
  
  # SUB: Y - D, on Y
  x n3
  x n4
  x n5
  x n6
  x n7
  cx d1,n4
  cx d2,n5
  cx d3,n6
  cx d4,n7
  cx d4,n8
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n3,d1
  toffoli d1,n4,d2
  toffoli d2,n5,d3
  toffoli d3,n6,d4
  toffoli d4,n7,n8
  cx d4,n7
  toffoli d3,n6,d4
  cx d3,n6
  toffoli d2,n5,d3
  cx d2,n5
  toffoli d1,n4,d2
  cx d1,n4
  toffoli d0,n3,d1
  cx d0,n3
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n7
  cx d3,n6
  cx d2,n5
  cx d1,n4
  x n3
  x n4
  x n5
  x n6
  x n7
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n4
  cx d2,n5
  cx d3,n6
  cx d4,n7
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n3,d1
  toffoli d1,n4,d2
  toffoli d2,n5,d3
  toffoli d3,n6,d4
  toffoli n8,d4,n7
  toffoli d3,n6,d4
  toffoli n8,d3,n6
  toffoli d2,n5,d3
  toffoli n8,d2,n5
  toffoli d1,n4,d2
  toffoli n8,d1,n4
  toffoli d0,n3,d1
  toffoli n8,d0,n3
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n7
  cx d3,n6
  cx d2,n5
  cx d1,n4
  
  # x on Z
  x n8
  
  ### Iteration 13, where Y = [n2, ..., n6] and Z = n7
  
  # SUB: Y - D, on Y
  x n2
  x n3
  x n4
  x n5
  x n6
  cx d1,n3
  cx d2,n4
  cx d3,n5
  cx d4,n6
  cx d4,n7
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n2,d1
  toffoli d1,n3,d2
  toffoli d2,n4,d3
  toffoli d3,n5,d4
  toffoli d4,n6,n7
  cx d4,n6
  toffoli d3,n5,d4
  cx d3,n5
  toffoli d2,n4,d3
  cx d2,n4
  toffoli d1,n3,d2
  cx d1,n3
  toffoli d0,n2,d1
  cx d0,n2
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n6
  cx d3,n5
  cx d2,n4
  cx d1,n3
  x n2
  x n3
  x n4
  x n5
  x n6
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n3
  cx d2,n4
  cx d3,n5
  cx d4,n6
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n2,d1
  toffoli d1,n3,d2
  toffoli d2,n4,d3
  toffoli d3,n5,d4
  toffoli n7,d4,n6
  toffoli d3,n5,d4
  toffoli n7,d3,n5
  toffoli d2,n4,d3
  toffoli n7,d2,n4
  toffoli d1,n3,d2
  toffoli n7,d1,n3
  toffoli d0,n2,d1
  toffoli n7,d0,n2
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n6
  cx d3,n5
  cx d2,n4
  cx d1,n3
  
  # x on Z
  x n7
  
  ### Iteration 14, where Y = [n1, ..., n5] and Z = n6
  
  # SUB: Y - D, on Y
  x n1
  x n2
  x n3
  x n4
  x n5
  cx d1,n2
  cx d2,n3
  cx d3,n4
  cx d4,n5
  cx d4,n6
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n1,d1
  toffoli d1,n2,d2
  toffoli d2,n3,d3
  toffoli d3,n4,d4
  toffoli d4,n5,n6
  cx d4,n5
  toffoli d3,n4,d4
  cx d3,n4
  toffoli d2,n3,d3
  cx d2,n3
  toffoli d1,n2,d2
  cx d1,n2
  toffoli d0,n1,d1
  cx d0,n1
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n5
  cx d3,n4
  cx d2,n3
  cx d1,n2
  x n1
  x n2
  x n3
  x n4
  x n5
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n2
  cx d2,n3
  cx d3,n4
  cx d4,n5
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n1,d1
  toffoli d1,n2,d2
  toffoli d2,n3,d3
  toffoli d3,n4,d4
  toffoli n6,d4,n5
  toffoli d3,n4,d4
  toffoli n6,d3,n4
  toffoli d2,n3,d3
  toffoli n6,d2,n3
  toffoli d1,n2,d2
  toffoli n6,d1,n2
  toffoli d0,n1,d1
  toffoli n6,d0,n1
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n5
  cx d3,n4
  cx d2,n3
  cx d1,n2
  
  # x on Z
  x n6
  
  ### Iteration 15, where Y = [n0, ..., n4] and Z = n5
  
  # SUB: Y - D, on Y
  x n0
  x n1
  x n2
  x n3
  x n4
  cx d1,n1
  cx d2,n2
  cx d3,n3
  cx d4,n4
  cx d4,n5
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n0,d1
  toffoli d1,n1,d2
  toffoli d2,n2,d3
  toffoli d3,n3,d4
  toffoli d4,n4,n5
  cx d4,n4
  toffoli d3,n3,d4
  cx d3,n3
  toffoli d2,n2,d3
  cx d2,n2
  toffoli d1,n1,d2
  cx d1,n1
  toffoli d0,n0,d1
  cx d0,n0
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n4
  cx d3,n3
  cx d2,n2
  cx d1,n1
  x n0
  x n1
  x n2
  x n3
  x n4
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx d1,n1
  cx d2,n2
  cx d3,n3
  cx d4,n4
  cx d3,d4
  cx d2,d3
  cx d1,d2
  toffoli d0,n0,d1
  toffoli d1,n1,d2
  toffoli d2,n2,d3
  toffoli d3,n3,d4
  toffoli n5,d4,n4
  toffoli d3,n3,d4
  toffoli n5,d3,n3
  toffoli d2,n2,d3
  toffoli n5,d2,n2
  toffoli d1,n1,d2
  toffoli n5,d1,n1
  toffoli d0,n0,d1
  toffoli n5,d0,n0
  cx d1,d2
  cx d2,d3
  cx d3,d4
  cx d4,n4
  cx d3,n3
  cx d2,n2
  cx d1,n1
  
  # x on Z
  x n5

.result
  measure
  display