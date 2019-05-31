#function: Ry(c*(x^k) rotation for c=0.25, k=3)

qubits 9

.init
  map q0,c0
  map q1,c1
  map q2,c2
  map q3,c3
  map q4,c4
  map q5,a0
  map q6,a1
  map q7,a2
  map q8,r
  x c0
  x c2
  x c4
  display

.ry_cx_to_the_k
  ### Rotations controlled by 0 qubit(s)
  
  # Rotation c*theta with c = 0.25 and theta = 1.0, controlled by [0]
  ry r,0.0078125
  cx c0,r
  ry r,-0.015625
  cx c0,r
  ry r,0.0078125
  
  # Rotation c*theta with c = 0.25 and theta = 0.125, controlled by [1]
  ry r,0.0009765625
  cx c1,r
  ry r,-0.001953125
  cx c1,r
  ry r,0.0009765625
  
  # Rotation c*theta with c = 0.25 and theta = 0.015625, controlled by [2]
  ry r,0.0001220703125
  cx c2,r
  ry r,-0.000244140625
  cx c2,r
  ry r,0.0001220703125
  
  # Rotation c*theta with c = 0.25 and theta = 0.001953125, controlled by [3]
  ry r,1.52587890625e-05
  cx c3,r
  ry r,-3.0517578125e-05
  cx c3,r
  ry r,1.52587890625e-05
  
  # Rotation c*theta with c = 0.25 and theta = 0.000244140625, controlled by [4]
  ry r,1.9073486328125e-06
  cx c4,r
  ry r,-3.814697265625e-06
  cx c4,r
  ry r,1.9073486328125e-06
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.25 and theta = 2.25, controlled by [0, 1]
  ry r,0.017578125
  cx c1,r
  ry r,-0.017578125
  cx c0,c1
  cx c0,r
  ry r,0.017578125
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.017578125
  cx c0,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.28125, controlled by [1, 2]
  ry r,0.002197265625
  cx c2,r
  ry r,-0.002197265625
  cx c1,c2
  cx c1,r
  ry r,0.002197265625
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.002197265625
  cx c1,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.1171875, controlled by [1, 3]
  ry r,0.00091552734375
  cx c3,r
  ry r,-0.00091552734375
  cx c1,c3
  cx c1,r
  ry r,0.00091552734375
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.00091552734375
  cx c1,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.052734375, controlled by [1, 4]
  ry r,0.0004119873046875
  cx c4,r
  ry r,-0.0004119873046875
  cx c1,c4
  cx c1,r
  ry r,0.0004119873046875
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.0004119873046875
  cx c1,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.0146484375, controlled by [2, 4]
  ry r,0.00011444091796875
  cx c4,r
  ry r,-0.00011444091796875
  cx c2,c4
  cx c2,r
  ry r,0.00011444091796875
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-0.00011444091796875
  cx c2,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.03515625, controlled by [2, 3]
  ry r,0.000274658203125
  cx c3,r
  ry r,-0.000274658203125
  cx c2,c3
  cx c2,r
  ry r,0.000274658203125
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-0.000274658203125
  cx c2,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.19921875, controlled by [0, 4]
  ry r,0.001556396484375
  cx c4,r
  ry r,-0.001556396484375
  cx c0,c4
  cx c0,r
  ry r,0.001556396484375
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.001556396484375
  cx c0,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.421875, controlled by [0, 3]
  ry r,0.0032958984375
  cx c3,r
  ry r,-0.0032958984375
  cx c0,c3
  cx c0,r
  ry r,0.0032958984375
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.0032958984375
  cx c0,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.00439453125, controlled by [3, 4]
  ry r,3.4332275390625e-05
  cx c4,r
  ry r,-3.4332275390625e-05
  cx c3,c4
  cx c3,r
  ry r,3.4332275390625e-05
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-3.4332275390625e-05
  cx c3,r
  
  # Rotation c*theta with c = 0.25 and theta = 0.9375, controlled by [0, 2]
  ry r,0.00732421875
  cx c2,r
  ry r,-0.00732421875
  cx c0,c2
  cx c0,r
  ry r,0.00732421875
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.00732421875
  cx c0,r
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.25 and theta = 0.046875, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.0003662109375
  cx a0,r
  ry r,-0.0003662109375
  cx c4,a0
  cx c4,r
  ry r,0.0003662109375
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0003662109375
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.75, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,0.005859375
  cx a0,r
  ry r,-0.005859375
  cx c2,a0
  cx c2,r
  ry r,0.005859375
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-0.005859375
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.1875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.00146484375
  cx a0,r
  ry r,-0.00146484375
  cx c3,a0
  cx c3,r
  ry r,0.00146484375
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.00146484375
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.0029296875
  cx a0,r
  ry r,-0.0029296875
  cx c3,a0
  cx c3,r
  ry r,0.0029296875
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.0029296875
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.09375, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.000732421875
  cx a0,r
  ry r,-0.000732421875
  cx c3,a0
  cx c3,r
  ry r,0.000732421875
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.000732421875
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.0234375, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.00018310546875
  cx a0,r
  ry r,-0.00018310546875
  cx c4,a0
  cx c4,r
  ry r,0.00018310546875
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.00018310546875
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.1875, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.00146484375
  cx a0,r
  ry r,-0.00146484375
  cx c4,a0
  cx c4,r
  ry r,0.00146484375
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.00146484375
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.046875, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.0003662109375
  cx a0,r
  ry r,-0.0003662109375
  cx c4,a0
  cx c4,r
  ry r,0.0003662109375
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0003662109375
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.01171875, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,9.1552734375e-05
  cx a0,r
  ry r,-9.1552734375e-05
  cx c4,a0
  cx c4,r
  ry r,9.1552734375e-05
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-9.1552734375e-05
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.25 and theta = 0.09375, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.000732421875
  cx a0,r
  ry r,-0.000732421875
  cx c4,a0
  cx c4,r
  ry r,0.000732421875
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.000732421875
  cx c4,r
  toffoli c0,c2,a0
  

.result
  display