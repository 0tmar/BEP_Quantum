#function: Ancilla rotation for c=1, approximated using the m=5 first terms of the taylor expansion of arcsin(x)

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
  # x c1
  # x c2
  # x c3
  # x c4
  display

.ry_c_x_to_the_k
  ##### Performs Ry(c*(x^k)) for c=1.0 and k=1
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 1.0 and theta = 1.0, controlled by [0]
  ry r,0.5
  cx c0,r
  ry r,-1.0
  cx c0,r
  ry r,0.5
  
  # Rotation c*theta with c = 1.0 and theta = 0.5, controlled by [1]
  ry r,0.25
  cx c1,r
  ry r,-0.5
  cx c1,r
  ry r,0.25
  
  # Rotation c*theta with c = 1.0 and theta = 0.25, controlled by [2]
  ry r,0.125
  cx c2,r
  ry r,-0.25
  cx c2,r
  ry r,0.125
  
  # Rotation c*theta with c = 1.0 and theta = 0.125, controlled by [3]
  ry r,0.0625
  cx c3,r
  ry r,-0.125
  cx c3,r
  ry r,0.0625
  
  # Rotation c*theta with c = 1.0 and theta = 0.0625, controlled by [4]
  ry r,0.03125
  cx c4,r
  ry r,-0.0625
  cx c4,r
  ry r,0.03125
  

.ry_c_x_to_the_k
  ##### Performs Ry(c*(x^k)) for c=0.16666666666666666 and k=3
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 1.0, controlled by [0]
  ry r,0.08333333333333333
  cx c0,r
  ry r,-0.16666666666666666
  cx c0,r
  ry r,0.08333333333333333
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.125, controlled by [1]
  ry r,0.010416666666666666
  cx c1,r
  ry r,-0.020833333333333332
  cx c1,r
  ry r,0.010416666666666666
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.015625, controlled by [2]
  ry r,0.0013020833333333333
  cx c2,r
  ry r,-0.0026041666666666665
  cx c2,r
  ry r,0.0013020833333333333
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.001953125, controlled by [3]
  ry r,0.00016276041666666666
  cx c3,r
  ry r,-0.0003255208333333333
  cx c3,r
  ry r,0.00016276041666666666
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.000244140625, controlled by [4]
  ry r,2.0345052083333332e-05
  cx c4,r
  ry r,-4.0690104166666664e-05
  cx c4,r
  ry r,2.0345052083333332e-05
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 2.25, controlled by [0, 1]
  ry r,0.1875
  cx c1,r
  ry r,-0.1875
  cx c0,c1
  cx c0,r
  ry r,0.1875
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.1875
  cx c0,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.28125, controlled by [1, 2]
  ry r,0.0234375
  cx c2,r
  ry r,-0.0234375
  cx c1,c2
  cx c1,r
  ry r,0.0234375
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.0234375
  cx c1,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.1171875, controlled by [1, 3]
  ry r,0.009765625
  cx c3,r
  ry r,-0.009765625
  cx c1,c3
  cx c1,r
  ry r,0.009765625
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.009765625
  cx c1,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.052734375, controlled by [1, 4]
  ry r,0.00439453125
  cx c4,r
  ry r,-0.00439453125
  cx c1,c4
  cx c1,r
  ry r,0.00439453125
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.00439453125
  cx c1,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.0146484375, controlled by [2, 4]
  ry r,0.001220703125
  cx c4,r
  ry r,-0.001220703125
  cx c2,c4
  cx c2,r
  ry r,0.001220703125
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-0.001220703125
  cx c2,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.03515625, controlled by [2, 3]
  ry r,0.0029296875
  cx c3,r
  ry r,-0.0029296875
  cx c2,c3
  cx c2,r
  ry r,0.0029296875
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-0.0029296875
  cx c2,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.19921875, controlled by [0, 4]
  ry r,0.0166015625
  cx c4,r
  ry r,-0.0166015625
  cx c0,c4
  cx c0,r
  ry r,0.0166015625
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.0166015625
  cx c0,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.421875, controlled by [0, 3]
  ry r,0.03515625
  cx c3,r
  ry r,-0.03515625
  cx c0,c3
  cx c0,r
  ry r,0.03515625
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.03515625
  cx c0,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.00439453125, controlled by [3, 4]
  ry r,0.0003662109375
  cx c4,r
  ry r,-0.0003662109375
  cx c3,c4
  cx c3,r
  ry r,0.0003662109375
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-0.0003662109375
  cx c3,r
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.9375, controlled by [0, 2]
  ry r,0.078125
  cx c2,r
  ry r,-0.078125
  cx c0,c2
  cx c0,r
  ry r,0.078125
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.078125
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.046875, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.00390625
  cx a0,r
  ry r,-0.00390625
  cx c4,a0
  cx c4,r
  ry r,0.00390625
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.00390625
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.75, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,0.0625
  cx a0,r
  ry r,-0.0625
  cx c2,a0
  cx c2,r
  ry r,0.0625
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-0.0625
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.1875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.015625
  cx a0,r
  ry r,-0.015625
  cx c3,a0
  cx c3,r
  ry r,0.015625
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.015625
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.03125
  cx a0,r
  ry r,-0.03125
  cx c3,a0
  cx c3,r
  ry r,0.03125
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.03125
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.09375, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.0078125
  cx a0,r
  ry r,-0.0078125
  cx c3,a0
  cx c3,r
  ry r,0.0078125
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.0078125
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.0234375, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.001953125
  cx a0,r
  ry r,-0.001953125
  cx c4,a0
  cx c4,r
  ry r,0.001953125
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.001953125
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.1875, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.015625
  cx a0,r
  ry r,-0.015625
  cx c4,a0
  cx c4,r
  ry r,0.015625
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.015625
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.046875, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.00390625
  cx a0,r
  ry r,-0.00390625
  cx c4,a0
  cx c4,r
  ry r,0.00390625
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.00390625
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.01171875, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,0.0009765625
  cx a0,r
  ry r,-0.0009765625
  cx c4,a0
  cx c4,r
  ry r,0.0009765625
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0009765625
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.16666666666666666 and theta = 0.09375, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.0078125
  cx a0,r
  ry r,-0.0078125
  cx c4,a0
  cx c4,r
  ry r,0.0078125
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0078125
  cx c4,r
  toffoli c0,c2,a0
  

.ry_c_x_to_the_k
  ##### Performs Ry(c*(x^k)) for c=0.075 and k=5
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.075 and theta = 1.0, controlled by [0]
  ry r,0.0375
  cx c0,r
  ry r,-0.075
  cx c0,r
  ry r,0.0375
  
  # Rotation c*theta with c = 0.075 and theta = 0.03125, controlled by [1]
  ry r,0.001171875
  cx c1,r
  ry r,-0.00234375
  cx c1,r
  ry r,0.001171875
  
  # Rotation c*theta with c = 0.075 and theta = 0.0009765625, controlled by [2]
  ry r,3.662109375e-05
  cx c2,r
  ry r,-7.32421875e-05
  cx c2,r
  ry r,3.662109375e-05
  
  # Rotation c*theta with c = 0.075 and theta = 3.0517578125e-05, controlled by [3]
  ry r,1.1444091796875e-06
  cx c3,r
  ry r,-2.288818359375e-06
  cx c3,r
  ry r,1.1444091796875e-06
  
  # Rotation c*theta with c = 0.075 and theta = 9.5367431640625e-07, controlled by [4]
  ry r,3.5762786865234374e-08
  cx c4,r
  ry r,-7.152557373046875e-08
  cx c4,r
  ry r,3.5762786865234374e-08
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.075 and theta = 6.5625, controlled by [0, 1]
  ry r,0.24609375
  cx c1,r
  ry r,-0.24609375
  cx c0,c1
  cx c0,r
  ry r,0.24609375
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.24609375
  cx c0,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.205078125, controlled by [1, 2]
  ry r,0.0076904296875
  cx c2,r
  ry r,-0.0076904296875
  cx c1,c2
  cx c1,r
  ry r,0.0076904296875
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.0076904296875
  cx c1,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.0640869140625, controlled by [1, 3]
  ry r,0.00240325927734375
  cx c3,r
  ry r,-0.00240325927734375
  cx c1,c3
  cx c1,r
  ry r,0.00240325927734375
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.00240325927734375
  cx c1,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.02506256103515625, controlled by [1, 4]
  ry r,0.0009398460388183594
  cx c4,r
  ry r,-0.0009398460388183594
  cx c1,c4
  cx c1,r
  ry r,0.0009398460388183594
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.0009398460388183594
  cx c1,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.002002716064453125, controlled by [2, 4]
  ry r,7.510185241699219e-05
  cx c4,r
  ry r,-7.510185241699219e-05
  cx c2,c4
  cx c2,r
  ry r,7.510185241699219e-05
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-7.510185241699219e-05
  cx c2,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.00640869140625, controlled by [2, 3]
  ry r,0.000240325927734375
  cx c3,r
  ry r,-0.000240325927734375
  cx c2,c3
  cx c2,r
  ry r,0.000240325927734375
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-0.000240325927734375
  cx c2,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.3540802001953125, controlled by [0, 4]
  ry r,0.013278007507324219
  cx c4,r
  ry r,-0.013278007507324219
  cx c0,c4
  cx c0,r
  ry r,0.013278007507324219
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.013278007507324219
  cx c0,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.802001953125, controlled by [0, 3]
  ry r,0.0300750732421875
  cx c3,r
  ry r,-0.0300750732421875
  cx c0,c3
  cx c0,r
  ry r,0.0300750732421875
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.0300750732421875
  cx c0,r
  
  # Rotation c*theta with c = 0.075 and theta = 0.0002002716064453125, controlled by [3, 4]
  ry r,7.510185241699219e-06
  cx c4,r
  ry r,-7.510185241699219e-06
  cx c3,c4
  cx c3,r
  ry r,7.510185241699219e-06
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-7.510185241699219e-06
  cx c3,r
  
  # Rotation c*theta with c = 0.075 and theta = 2.05078125, controlled by [0, 2]
  ry r,0.076904296875
  cx c2,r
  ry r,-0.076904296875
  cx c0,c2
  cx c0,r
  ry r,0.076904296875
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.076904296875
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.075 and theta = 0.205078125, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.0076904296875
  cx a0,r
  ry r,-0.0076904296875
  cx c4,a0
  cx c4,r
  ry r,0.0076904296875
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0076904296875
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.075 and theta = 6.5625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,0.24609375
  cx a0,r
  ry r,-0.24609375
  cx c2,a0
  cx c2,r
  ry r,0.24609375
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-0.24609375
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.075 and theta = 1.0546875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.03955078125
  cx a0,r
  ry r,-0.03955078125
  cx c3,a0
  cx c3,r
  ry r,0.03955078125
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.03955078125
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.075 and theta = 2.87109375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.107666015625
  cx a0,r
  ry r,-0.107666015625
  cx c3,a0
  cx c3,r
  ry r,0.107666015625
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.107666015625
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.205078125, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.0076904296875
  cx a0,r
  ry r,-0.0076904296875
  cx c3,a0
  cx c3,r
  ry r,0.0076904296875
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.0076904296875
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.032958984375, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.0012359619140625
  cx a0,r
  ry r,-0.0012359619140625
  cx c4,a0
  cx c4,r
  ry r,0.0012359619140625
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0012359619140625
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.075 and theta = 1.34033203125, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.050262451171875
  cx a0,r
  ry r,-0.050262451171875
  cx c4,a0
  cx c4,r
  ry r,0.050262451171875
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.050262451171875
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.0897216796875, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.00336456298828125
  cx a0,r
  ry r,-0.00336456298828125
  cx c4,a0
  cx c4,r
  ry r,0.00336456298828125
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.00336456298828125
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.00640869140625, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,0.000240325927734375
  cx a0,r
  ry r,-0.000240325927734375
  cx c4,a0
  cx c4,r
  ry r,0.000240325927734375
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.000240325927734375
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.487060546875, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.0182647705078125
  cx a0,r
  ry r,-0.0182647705078125
  cx c4,a0
  cx c4,r
  ry r,0.0182647705078125
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0182647705078125
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.075 and theta = 0.054931640625, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.0020599365234375
  cx a1,r
  ry r,-0.0020599365234375
  cx c4,a1
  cx c4,r
  ry r,0.0020599365234375
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.0020599365234375
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.3955078125, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,0.01483154296875
  cx a1,r
  ry r,-0.01483154296875
  cx c4,a1
  cx c4,r
  ry r,0.01483154296875
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.01483154296875
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.849609375, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.0318603515625
  cx a1,r
  ry r,-0.0318603515625
  cx c4,a1
  cx c4,r
  ry r,0.0318603515625
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.0318603515625
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.075 and theta = 0.16845703125, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,0.006317138671875
  cx a1,r
  ry r,-0.006317138671875
  cx c4,a1
  cx c4,r
  ry r,0.006317138671875
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.006317138671875
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.075 and theta = 1.7578125, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.06591796875
  cx a1,r
  ry r,-0.06591796875
  cx c3,a1
  cx c3,r
  ry r,0.06591796875
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-0.06591796875
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 0.075 and theta = 0.1171875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,0.00439453125
  cx a2,r
  ry r,-0.00439453125
  cx c4,a2
  cx c4,r
  ry r,0.00439453125
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-0.00439453125
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.ry_c_x_to_the_k
  ##### Performs Ry(c*(x^k)) for c=0.044642857142857144 and k=7
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 1.0, controlled by [0]
  ry r,0.022321428571428572
  cx c0,r
  ry r,-0.044642857142857144
  cx c0,r
  ry r,0.022321428571428572
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.0078125, controlled by [1]
  ry r,0.00017438616071428572
  cx c1,r
  ry r,-0.00034877232142857144
  cx c1,r
  ry r,0.00017438616071428572
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 6.103515625e-05, controlled by [2]
  ry r,1.3623918805803572e-06
  cx c2,r
  ry r,-2.7247837611607143e-06
  cx c2,r
  ry r,1.3623918805803572e-06
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 4.76837158203125e-07, controlled by [3]
  ry r,1.064368656703404e-08
  cx c3,r
  ry r,-2.128737313406808e-08
  cx c3,r
  ry r,1.064368656703404e-08
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 3.725290298461914e-09, controlled by [4]
  ry r,8.315380130495344e-11
  cx c4,r
  ry r,-1.6630760260990688e-10
  cx c4,r
  ry r,8.315380130495344e-11
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 16.078125, controlled by [0, 1]
  ry r,0.35888671875
  cx c1,r
  ry r,-0.35888671875
  cx c0,c1
  cx c0,r
  ry r,0.35888671875
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.35888671875
  cx c0,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.1256103515625, controlled by [1, 2]
  ry r,0.002803802490234375
  cx c2,r
  ry r,-0.002803802490234375
  cx c1,c2
  cx c1,r
  ry r,0.002803802490234375
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.002803802490234375
  cx c1,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.029439926147460938, controlled by [1, 3]
  ry r,0.0006571412086486816
  cx c3,r
  ry r,-0.0006571412086486816
  cx c1,c3
  cx c1,r
  ry r,0.0006571412086486816
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.0006571412086486816
  cx c1,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.010005444288253784, controlled by [1, 4]
  ry r,0.00022333581000566483
  cx c4,r
  ry r,-0.00022333581000566483
  cx c1,c4
  cx c1,r
  ry r,0.00022333581000566483
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.00022333581000566483
  cx c1,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.00022999942302703857, controlled by [2, 4]
  ry r,5.133915692567825e-06
  cx c4,r
  ry r,-5.133915692567825e-06
  cx c2,c4
  cx c2,r
  ry r,5.133915692567825e-06
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-5.133915692567825e-06
  cx c2,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.0009813308715820312, controlled by [2, 3]
  ry r,2.1904706954956055e-05
  cx c3,r
  ry r,-2.1904706954956055e-05
  cx c2,c3
  cx c2,r
  ry r,2.1904706954956055e-05
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-2.1904706954956055e-05
  cx c2,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.5286306738853455, controlled by [0, 4]
  ry r,0.01179979182779789
  cx c4,r
  ry r,-0.01179979182779789
  cx c0,c4
  cx c0,r
  ry r,0.01179979182779789
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.01179979182779789
  cx c0,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 1.2806968688964844, controlled by [0, 3]
  ry r,0.028586983680725098
  cx c3,r
  ry r,-0.028586983680725098
  cx c0,c3
  cx c0,r
  ry r,0.028586983680725098
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.028586983680725098
  cx c0,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 7.666647434234619e-06, controlled by [3, 4]
  ry r,1.7113052308559418e-07
  cx c4,r
  ry r,-1.7113052308559418e-07
  cx c3,c4
  cx c3,r
  ry r,1.7113052308559418e-07
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-1.7113052308559418e-07
  cx c3,r
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 3.768310546875, controlled by [0, 2]
  ry r,0.08411407470703125
  cx c2,r
  ry r,-0.08411407470703125
  cx c0,c2
  cx c0,r
  ry r,0.08411407470703125
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.08411407470703125
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.5205960273742676, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.011620447039604187
  cx a0,r
  ry r,-0.011620447039604187
  cx c4,a0
  cx c4,r
  ry r,0.011620447039604187
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.011620447039604187
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 29.28515625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,0.6536865234375
  cx a0,r
  ry r,-0.6536865234375
  cx c2,a0
  cx c2,r
  ry r,0.6536865234375
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-0.6536865234375
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 3.242156982421875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.07236957550048828
  cx a0,r
  ry r,-0.07236957550048828
  cx c3,a0
  cx c3,r
  ry r,0.07236957550048828
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.07236957550048828
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 11.524749755859375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.2572488784790039
  cx a0,r
  ry r,-0.2572488784790039
  cx c3,a0
  cx c3,r
  ry r,0.2572488784790039
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.2572488784790039
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.228790283203125, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.005106925964355469
  cx a0,r
  ry r,-0.005106925964355469
  cx c3,a0
  cx c3,r
  ry r,0.005106925964355469
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.005106925964355469
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.0253293514251709, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.0005653873085975647
  cx a0,r
  ry r,-0.0005653873085975647
  cx c4,a0
  cx c4,r
  ry r,0.0005653873085975647
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0005653873085975647
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 5.112793922424316, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.11412486433982849
  cx a0,r
  ry r,-0.11412486433982849
  cx c4,a0
  cx c4,r
  ry r,0.11412486433982849
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.11412486433982849
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.09003710746765137, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.002009756863117218
  cx a0,r
  ry r,-0.002009756863117218
  cx c4,a0
  cx c4,r
  ry r,0.002009756863117218
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.002009756863117218
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.001787424087524414, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,3.98978590965271e-05
  cx a0,r
  ry r,-3.98978590965271e-05
  cx c4,a0
  cx c4,r
  ry r,3.98978590965271e-05
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-3.98978590965271e-05
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 1.4123454093933105, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.031525567173957825
  cx a0,r
  ry r,-0.031525567173957825
  cx c4,a0
  cx c4,r
  ry r,0.031525567173957825
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.031525567173957825
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.11640787124633789, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.002598389983177185
  cx a1,r
  ry r,-0.002598389983177185
  cx c4,a1
  cx c4,r
  ry r,0.002598389983177185
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.002598389983177185
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 2.8496646881103516, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,0.06360858678817749
  cx a1,r
  ry r,-0.06360858678817749
  cx c4,a1
  cx c4,r
  ry r,0.06360858678817749
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.06360858678817749
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 6.841678619384766, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.1527160406112671
  cx a1,r
  ry r,-0.1527160406112671
  cx c4,a1
  cx c4,r
  ry r,0.1527160406112671
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.1527160406112671
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 0.9281587600708008, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,0.02071782946586609
  cx a1,r
  ry r,-0.02071782946586609
  cx c4,a1
  cx c4,r
  ry r,0.02071782946586609
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.02071782946586609
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 14.90020751953125, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.3325939178466797
  cx a1,r
  ry r,-0.3325939178466797
  cx c3,a1
  cx c3,r
  ry r,0.3325939178466797
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-0.3325939178466797
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 0.044642857142857144 and theta = 2.58270263671875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,0.05764961242675781
  cx a2,r
  ry r,-0.05764961242675781
  cx c4,a2
  cx c4,r
  ry r,0.05764961242675781
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-0.05764961242675781
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.ry_c_x_to_the_k
  ##### Performs Ry(c*(x^k)) for c=0.030381944444444444 and k=9
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 1.0, controlled by [0]
  ry r,0.015190972222222222
  cx c0,r
  ry r,-0.030381944444444444
  cx c0,r
  ry r,0.015190972222222222
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.001953125, controlled by [1]
  ry r,2.9669867621527777e-05
  cx c1,r
  ry r,-5.9339735243055555e-05
  cx c1,r
  ry r,2.9669867621527777e-05
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 3.814697265625e-06, controlled by [2]
  ry r,5.794896019829644e-08
  cx c2,r
  ry r,-1.1589792039659288e-07
  cx c2,r
  ry r,5.794896019829644e-08
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 7.450580596923828e-09, controlled by [3]
  ry r,1.1318156288729773e-10
  cx c3,r
  ry r,-2.2636312577459547e-10
  cx c3,r
  ry r,1.1318156288729773e-10
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 1.4551915228366852e-11, controlled by [4]
  ry r,2.210577400142534e-13
  cx c4,r
  ry r,-4.421154800285068e-13
  cx c4,r
  ry r,2.210577400142534e-13
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 37.44140625, controlled by [0, 1]
  ry r,0.5687713623046875
  cx c1,r
  ry r,-0.5687713623046875
  cx c0,c1
  cx c0,r
  ry r,0.5687713623046875
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.5687713623046875
  cx c0,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.07312774658203125, controlled by [1, 2]
  ry r,0.0011108815670013428
  cx c2,r
  ry r,-0.0011108815670013428
  cx c1,c2
  cx c1,r
  ry r,0.0011108815670013428
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.0011108815670013428
  cx c1,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.012598782777786255, controlled by [1, 3]
  ry r,0.00019138775921116272
  cx c3,r
  ry r,-0.00019138775921116272
  cx c1,c3
  cx c1,r
  ry r,0.00019138775921116272
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.00019138775921116272
  cx c1,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.003684585099108517, controlled by [1, 4]
  ry r,5.59724298909714e-05
  cx c4,r
  ry r,-5.59724298909714e-05
  cx c1,c4
  cx c1,r
  ry r,5.59724298909714e-05
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-5.59724298909714e-05
  cx c1,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 2.460699761286378e-05, controlled by [2, 4]
  ry r,3.738042172093022e-07
  cx c4,r
  ry r,-3.738042172093022e-07
  cx c2,c4
  cx c2,r
  ry r,3.738042172093022e-07
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-3.738042172093022e-07
  cx c2,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.00014282763004302979, controlled by [2, 3]
  ry r,2.1696905605494976e-06
  cx c3,r
  ry r,-2.1696905605494976e-06
  cx c2,c3
  cx c2,r
  ry r,2.1696905605494976e-06
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-2.1696905605494976e-06
  cx c2,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.7256807258818299, controlled by [0, 4]
  ry r,0.011023795749072937
  cx c4,r
  ry r,-0.011023795749072937
  cx c0,c4
  cx c0,r
  ry r,0.011023795749072937
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.011023795749072937
  cx c0,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 1.8865075707435608, controlled by [0, 3]
  ry r,0.028657884104177356
  cx c3,r
  ry r,-0.028657884104177356
  cx c0,c3
  cx c0,r
  ry r,0.028657884104177356
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.028657884104177356
  cx c0,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 2.7896021492779255e-07, controlled by [3, 4]
  ry r,4.2376768760732375e-09
  cx c4,r
  ry r,-4.2376768760732375e-09
  cx c3,c4
  cx c3,r
  ry r,4.2376768760732375e-09
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-4.2376768760732375e-09
  cx c3,r
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 6.4505767822265625, controlled by [0, 2]
  ry r,0.09799053271611531
  cx c2,r
  ry r,-0.09799053271611531
  cx c0,c2
  cx c0,r
  ry r,0.09799053271611531
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.09799053271611531
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 1.0835354402661324, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.016459956774876144
  cx a0,r
  ry r,-0.016459956774876144
  cx c4,a0
  cx c4,r
  ry r,0.016459956774876144
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.016459956774876144
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 108.9697265625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,1.6553560892740886
  cx a0,r
  ry r,-1.6553560892740886
  cx c2,a0
  cx c2,r
  ry r,1.6553560892740886
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-1.6553560892740886
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 8.23084831237793, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.12503458807865778
  cx a0,r
  ry r,-0.12503458807865778
  cx c3,a0
  cx c3,r
  ry r,0.12503458807865778
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.12503458807865778
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 38.66721153259277, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.5873925363024076
  cx a0,r
  ry r,-0.5873925363024076
  cx c3,a0
  cx c3,r
  ry r,0.5873925363024076
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.5873925363024076
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.2128314971923828, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.0032331173618634543
  cx a0,r
  ry r,-0.0032331173618634543
  cx c3,a0
  cx c3,r
  ry r,0.0032331173618634543
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.0032331173618634543
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.016075875610113144, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.00024420817984112847
  cx a0,r
  ry r,-0.00024420817984112847
  cx c4,a0
  cx c4,r
  ry r,0.00024420817984112847
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.00024420817984112847
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 16.338426545262337, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.24819658380389834
  cx a0,r
  ry r,-0.24819658380389834
  cx c4,a0
  cx c4,r
  ry r,0.24819658380389834
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.24819658380389834
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.07552189752459526, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.0011472510474656399
  cx a0,r
  ry r,-0.0011472510474656399
  cx c4,a0
  cx c4,r
  ry r,0.0011472510474656399
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0011472510474656399
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.0004156865179538727, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,6.314682347389559e-06
  cx a0,r
  ry r,-6.314682347389559e-06
  cx c4,a0
  cx c4,r
  ry r,6.314682347389559e-06
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-6.314682347389559e-06
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 3.3820099756121635, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.051376019594802834
  cx a0,r
  ry r,-0.051376019594802834
  cx c4,a0
  cx c4,r
  ry r,0.051376019594802834
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.051376019594802834
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 0.163043774664402, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.0024767934519331902
  cx a1,r
  ry r,-0.0024767934519331902
  cx c4,a1
  cx c4,r
  ry r,0.0024767934519331902
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.0024767934519331902
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 13.789967447519302, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,0.2094830124406144
  cx a1,r
  ry r,-0.2094830124406144
  cx c4,a1
  cx c4,r
  ry r,0.2094830124406144
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.2094830124406144
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 36.644615828990936, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.5566673411522061
  cx a1,r
  ry r,-0.5566673411522061
  cx c4,a1
  cx c4,r
  ry r,0.5566673411522061
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.5566673411522061
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 3.450474813580513, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,0.052416067046578974
  cx a1,r
  ry r,-0.052416067046578974
  cx c4,a1
  cx c4,r
  ry r,0.052416067046578974
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.052416067046578974
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 83.47841262817383, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,1.2681182473897934
  cx a1,r
  ry r,-1.2681182473897934
  cx c3,a1
  cx c3,r
  ry r,1.2681182473897934
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-1.2681182473897934
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 0.030381944444444444 and theta = 22.6483154296875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,0.34404993057250977
  cx a2,r
  ry r,-0.34404993057250977
  cx c4,a2
  cx c4,r
  ry r,0.34404993057250977
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-0.34404993057250977
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.result
  display