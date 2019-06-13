#function: Complete QLSA for the Cao matrix

qubits 20

.init
  map q0,f
  map q1,e0
  map q2,e1
  map q3,d0
  map q4,d1
  map q5,d2
  map q6,d3
  map q7,d4
  map q8,c0
  map q9,c1
  map q10,c2
  map q11,c3
  map q12,c4
  map q13,b0
  map q14,b1
  map q15,b2
  map q16,b3
  map q17,b4
  map q18,a0
  map q19,a1
  h b0
  h b1
  h b2
  h b3
  h a0
  h a1
  x d4
  display

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/16) on a0 and a1, controlled by b0
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b0,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b0,a0
  tdag b0
  tdag a1
  cx b0,a1
  t a1
  cx b0,a1
  
  # cRx(a,c,-0.196)
  rx a1,-0.049
  cz b0,a1
  rx a1,0.098
  cz b0,a1
  rx a1,-0.049
  
  # cVdag(a,c)
  t b0
  x b0
  tdag b0
  x b0
  h a1
  tdag a1
  tdag a1
  cx b0,a1
  h a1
  t a1
  cx b0,a1
  t a1
  h a1
  s a1
  x a1
  tdag b0
  tdag a1
  cx b0,a1
  t a1
  cx b0,a1
  t b0
  x b0
  tdag b0
  x b0
  h a1
  tdag a1
  tdag a1
  cx b0,a1
  h a1
  t a1
  cx b0,a1
  t a1
  h a1
  s a1
  x a1
  
  # Rz(a,0.375)
  rz b0,0.375
  
  # cRx(a,b,-0.982)
  rx a0,-0.2455
  cz b0,a0
  rx a0,0.491
  cz b0,a0
  rx a0,-0.2455
  
  # Rz(a,1.883)
  rz b0,1.883
  
  # Toffoli(a,b,c)
  toffoli b0,a0,a1
  
  # cRx(a,b,-0.589)
  rx a0,-0.14725
  cz b0,a0
  rx a0,0.2945
  cz b0,a0
  rx a0,-0.14725
  
  # Toffoli(a,b,c)
  toffoli b0,a0,a1
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b0,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b0,a0
  tdag b0
  tdag a1
  cx b0,a1
  t a1
  cx b0,a1

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/8) on a0 and a1, controlled by b1
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b1,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b1,a0
  tdag b1
  tdag a1
  cx b1,a1
  t a1
  cx b1,a1
  
  # cRx(a,c,-0.392)
  rx a1,-0.098
  cz b1,a1
  rx a1,0.196
  cz b1,a1
  rx a1,-0.098
  
  # cX(a,c)
  cx b1,a1
  
  # Rz(a,0.75)
  rz b1,0.75
  
  # cRx(a,b,-1.964)
  rx a0,-0.491
  cz b1,a0
  rx a0,0.982
  cz b1,a0
  rx a0,-0.491
  
  # Rz(a,3.766)
  rz b1,3.766
  
  # Toffoli(a,b,c)
  toffoli b1,a0,a1
  
  # cRx(a,b,-1.178)
  rx a0,-0.2945
  cz b1,a0
  rx a0,0.589
  cz b1,a0
  rx a0,-0.2945
  
  # Toffoli(a,b,c)
  toffoli b1,a0,a1
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b1,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b1,a0
  tdag b1
  tdag a1
  cx b1,a1
  t a1
  cx b1,a1

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/4) on a0 and a1, controlled by b2
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b2,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b2,a0
  tdag b2
  tdag a1
  cx b2,a1
  t a1
  cx b2,a1
  
  # cRx(a,c,-0.784)
  rx a1,-0.196
  cz b2,a1
  rx a1,0.392
  cz b2,a1
  rx a1,-0.196
  
  # Rz(a,1.5)
  rz b2,1.5
  
  # cRx(a,b,-3.928)
  rx a0,-0.982
  cz b2,a0
  rx a0,1.964
  cz b2,a0
  rx a0,-0.982
  
  # Rz(a,7.532)
  rz b2,7.532
  
  # Toffoli(a,b,c)
  toffoli b2,a0,a1
  
  # cRx(a,b,-2.356)
  rx a0,-0.589
  cz b2,a0
  rx a0,1.178
  cz b2,a0
  rx a0,-0.589
  
  # Toffoli(a,b,c)
  toffoli b2,a0,a1
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b2,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b2,a0
  tdag b2
  tdag a1
  cx b2,a1
  t a1
  cx b2,a1

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/2) on a0 and a1, controlled by b3
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b3,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b3,a0
  tdag b3
  tdag a1
  cx b3,a1
  t a1
  cx b3,a1
  
  # cRx(a,c,-1.568)
  rx a1,-0.392
  cz b3,a1
  rx a1,0.784
  cz b3,a1
  rx a1,-0.392
  
  # Rz(a,3.0)
  rz b3,3.0
  
  # cRx(a,b,-7.856)
  rx a0,-1.964
  cz b3,a0
  rx a0,3.928
  cz b3,a0
  rx a0,-1.964
  
  # Rz(a,15.064)
  rz b3,15.064
  
  # Toffoli(a,b,c)
  toffoli b3,a0,a1
  
  # cRx(a,b,-4.712)
  rx a0,-1.178
  cz b3,a0
  rx a0,2.356
  cz b3,a0
  rx a0,-1.178
  
  # Toffoli(a,b,c)
  toffoli b3,a0,a1
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b3,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b3,a0
  tdag b3
  tdag a1
  cx b3,a1
  t a1
  cx b3,a1

.iqft
  h b3
  cr b3,b2
  cz b3,b2
  h b2
  cr b3,b1
  swap b3,b2
  cr b2,b1
  swap b3,b2
  cz b3,b1
  cr b2,b1
  cz b2,b1
  h b1
  cr b3,b0
  swap b3,b2
  cr b2,b0
  swap b3,b2
  swap b3,b1
  cr b1,b0
  swap b3,b1
  cz b3,b0
  cr b2,b0
  swap b2,b1
  cr b1,b0
  swap b2,b1
  cz b2,b0
  cr b1,b0
  cz b1,b0
  h b0

.reverse
  swap b0,b3
  swap b1,b2

.display
  display

.thapliyal_division
  ##### Thapliyal division with N = [d0, ..., d4], O = [c0, ..., c4] and D = [b0, ..., b4]
  
  ### Iteration 1, where Y = [d4, ..., c3] and Z = c4
  
  # SUB: Y - D, on Y
  x d4
  x c0
  x c1
  x c2
  x c3
  cx b1,c0
  cx b2,c1
  cx b3,c2
  cx b4,c3
  cx b4,c4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d4,b1
  toffoli b1,c0,b2
  toffoli b2,c1,b3
  toffoli b3,c2,b4
  toffoli b4,c3,c4
  cx b4,c3
  toffoli b3,c2,b4
  cx b3,c2
  toffoli b2,c1,b3
  cx b2,c1
  toffoli b1,c0,b2
  cx b1,c0
  toffoli b0,d4,b1
  cx b0,d4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c3
  cx b3,c2
  cx b2,c1
  cx b1,c0
  x d4
  x c0
  x c1
  x c2
  x c3
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,c0
  cx b2,c1
  cx b3,c2
  cx b4,c3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d4,b1
  toffoli b1,c0,b2
  toffoli b2,c1,b3
  toffoli b3,c2,b4
  toffoli c4,b4,c3
  toffoli b3,c2,b4
  toffoli c4,b3,c2
  toffoli b2,c1,b3
  toffoli c4,b2,c1
  toffoli b1,c0,b2
  toffoli c4,b1,c0
  toffoli b0,d4,b1
  toffoli c4,b0,d4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c3
  cx b3,c2
  cx b2,c1
  cx b1,c0
  
  # x on Z
  x c4
  
  ### Iteration 2, where Y = [d3, ..., c2] and Z = c3
  
  # SUB: Y - D, on Y
  x d3
  x d4
  x c0
  x c1
  x c2
  cx b1,d4
  cx b2,c0
  cx b3,c1
  cx b4,c2
  cx b4,c3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d3,b1
  toffoli b1,d4,b2
  toffoli b2,c0,b3
  toffoli b3,c1,b4
  toffoli b4,c2,c3
  cx b4,c2
  toffoli b3,c1,b4
  cx b3,c1
  toffoli b2,c0,b3
  cx b2,c0
  toffoli b1,d4,b2
  cx b1,d4
  toffoli b0,d3,b1
  cx b0,d3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c2
  cx b3,c1
  cx b2,c0
  cx b1,d4
  x d3
  x d4
  x c0
  x c1
  x c2
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,d4
  cx b2,c0
  cx b3,c1
  cx b4,c2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d3,b1
  toffoli b1,d4,b2
  toffoli b2,c0,b3
  toffoli b3,c1,b4
  toffoli c3,b4,c2
  toffoli b3,c1,b4
  toffoli c3,b3,c1
  toffoli b2,c0,b3
  toffoli c3,b2,c0
  toffoli b1,d4,b2
  toffoli c3,b1,d4
  toffoli b0,d3,b1
  toffoli c3,b0,d3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c2
  cx b3,c1
  cx b2,c0
  cx b1,d4
  
  # x on Z
  x c3
  
  ### Iteration 3, where Y = [d2, ..., c1] and Z = c2
  
  # SUB: Y - D, on Y
  x d2
  x d3
  x d4
  x c0
  x c1
  cx b1,d3
  cx b2,d4
  cx b3,c0
  cx b4,c1
  cx b4,c2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d2,b1
  toffoli b1,d3,b2
  toffoli b2,d4,b3
  toffoli b3,c0,b4
  toffoli b4,c1,c2
  cx b4,c1
  toffoli b3,c0,b4
  cx b3,c0
  toffoli b2,d4,b3
  cx b2,d4
  toffoli b1,d3,b2
  cx b1,d3
  toffoli b0,d2,b1
  cx b0,d2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c1
  cx b3,c0
  cx b2,d4
  cx b1,d3
  x d2
  x d3
  x d4
  x c0
  x c1
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,d3
  cx b2,d4
  cx b3,c0
  cx b4,c1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d2,b1
  toffoli b1,d3,b2
  toffoli b2,d4,b3
  toffoli b3,c0,b4
  toffoli c2,b4,c1
  toffoli b3,c0,b4
  toffoli c2,b3,c0
  toffoli b2,d4,b3
  toffoli c2,b2,d4
  toffoli b1,d3,b2
  toffoli c2,b1,d3
  toffoli b0,d2,b1
  toffoli c2,b0,d2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c1
  cx b3,c0
  cx b2,d4
  cx b1,d3
  
  # x on Z
  x c2
  
  ### Iteration 4, where Y = [d1, ..., c0] and Z = c1
  
  # SUB: Y - D, on Y
  x d1
  x d2
  x d3
  x d4
  x c0
  cx b1,d2
  cx b2,d3
  cx b3,d4
  cx b4,c0
  cx b4,c1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d1,b1
  toffoli b1,d2,b2
  toffoli b2,d3,b3
  toffoli b3,d4,b4
  toffoli b4,c0,c1
  cx b4,c0
  toffoli b3,d4,b4
  cx b3,d4
  toffoli b2,d3,b3
  cx b2,d3
  toffoli b1,d2,b2
  cx b1,d2
  toffoli b0,d1,b1
  cx b0,d1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c0
  cx b3,d4
  cx b2,d3
  cx b1,d2
  x d1
  x d2
  x d3
  x d4
  x c0
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,d2
  cx b2,d3
  cx b3,d4
  cx b4,c0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d1,b1
  toffoli b1,d2,b2
  toffoli b2,d3,b3
  toffoli b3,d4,b4
  toffoli c1,b4,c0
  toffoli b3,d4,b4
  toffoli c1,b3,d4
  toffoli b2,d3,b3
  toffoli c1,b2,d3
  toffoli b1,d2,b2
  toffoli c1,b1,d2
  toffoli b0,d1,b1
  toffoli c1,b0,d1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c0
  cx b3,d4
  cx b2,d3
  cx b1,d2
  
  # x on Z
  x c1
  
  ### Iteration 5, where Y = [d0, ..., d4] and Z = c0
  
  # SUB: Y - D, on Y
  x d0
  x d1
  x d2
  x d3
  x d4
  cx b1,d1
  cx b2,d2
  cx b3,d3
  cx b4,d4
  cx b4,c0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d0,b1
  toffoli b1,d1,b2
  toffoli b2,d2,b3
  toffoli b3,d3,b4
  toffoli b4,d4,c0
  cx b4,d4
  toffoli b3,d3,b4
  cx b3,d3
  toffoli b2,d2,b3
  cx b2,d2
  toffoli b1,d1,b2
  cx b1,d1
  toffoli b0,d0,b1
  cx b0,d0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d4
  cx b3,d3
  cx b2,d2
  cx b1,d1
  x d0
  x d1
  x d2
  x d3
  x d4
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,d1
  cx b2,d2
  cx b3,d3
  cx b4,d4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d0,b1
  toffoli b1,d1,b2
  toffoli b2,d2,b3
  toffoli b3,d3,b4
  toffoli c0,b4,d4
  toffoli b3,d3,b4
  toffoli c0,b3,d3
  toffoli b2,d2,b3
  toffoli c0,b2,d2
  toffoli b1,d1,b2
  toffoli c0,b1,d1
  toffoli b0,d0,b1
  toffoli c0,b0,d0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d4
  cx b3,d3
  cx b2,d2
  cx b1,d1
  
  # x on Z
  x c0

.display
  display

.ancilla_rotation
  ####### Ancilla rotation using the first m terms of arcsin(x) and c=0.7853981633974483
  
  ##### Performs Ry(c*(x^k)) for c=0.7853981633974483 and k=1
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.25, controlled by [2]
  ry f,0.09817477042468103
  cx c2,f
  ry f,-0.19634954084936207
  cx c2,f
  ry f,0.09817477042468103
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 1.0, controlled by [0]
  ry f,0.39269908169872414
  cx c4,f
  ry f,-0.7853981633974483
  cx c4,f
  ry f,0.39269908169872414
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.125, controlled by [3]
  ry f,0.04908738521234052
  cx c1,f
  ry f,-0.09817477042468103
  cx c1,f
  ry f,0.04908738521234052
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.5, controlled by [1]
  ry f,0.19634954084936207
  cx c3,f
  ry f,-0.39269908169872414
  cx c3,f
  ry f,0.19634954084936207
  
  ##### Performs Ry(c*(x^k)) for c=0.08074551218828077 and k=3
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.015625, controlled by [2]
  ry f,0.0006308243139709435
  cx c2,f
  ry f,-0.001261648627941887
  cx c2,f
  ry f,0.0006308243139709435
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 1.0, controlled by [0]
  ry f,0.040372756094140386
  cx c4,f
  ry f,-0.08074551218828077
  cx c4,f
  ry f,0.040372756094140386
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.001953125, controlled by [3]
  ry f,7.885303924636794e-05
  cx c1,f
  ry f,-0.00015770607849273588
  cx c1,f
  ry f,7.885303924636794e-05
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.125, controlled by [1]
  ry f,0.005046594511767548
  cx c3,f
  ry f,-0.010093189023535096
  cx c3,f
  ry f,0.005046594511767548
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 2.25, controlled by [0, 1]
  ry f,0.09083870121181586
  cx c3,f
  ry f,-0.09083870121181586
  cx c4,c3
  cx c4,f
  ry f,0.09083870121181586
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.09083870121181586
  cx c4,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.28125, controlled by [1, 2]
  ry f,0.011354837651476983
  cx c2,f
  ry f,-0.011354837651476983
  cx c3,c2
  cx c3,f
  ry f,0.011354837651476983
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-0.011354837651476983
  cx c3,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.1171875, controlled by [1, 3]
  ry f,0.0047311823547820766
  cx c1,f
  ry f,-0.0047311823547820766
  cx c3,c1
  cx c3,f
  ry f,0.0047311823547820766
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-0.0047311823547820766
  cx c3,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.03515625, controlled by [2, 3]
  ry f,0.0014193547064346229
  cx c1,f
  ry f,-0.0014193547064346229
  cx c2,c1
  cx c2,f
  ry f,0.0014193547064346229
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-0.0014193547064346229
  cx c2,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.421875, controlled by [0, 3]
  ry f,0.017032256477215475
  cx c1,f
  ry f,-0.017032256477215475
  cx c4,c1
  cx c4,f
  ry f,0.017032256477215475
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.017032256477215475
  cx c4,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.9375, controlled by [0, 2]
  ry f,0.03784945883825661
  cx c2,f
  ry f,-0.03784945883825661
  cx c4,c2
  cx c4,f
  ry f,0.03784945883825661
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.03784945883825661
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.75, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.03027956707060529
  cx e0,f
  ry f,-0.03027956707060529
  cx c2,e0
  cx c2,f
  ry f,0.03027956707060529
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.03027956707060529
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.1875, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.007569891767651322
  cx e0,f
  ry f,-0.007569891767651322
  cx c1,e0
  cx c1,f
  ry f,0.007569891767651322
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.007569891767651322
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.09375, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,0.003784945883825661
  cx e0,f
  ry f,-0.003784945883825661
  cx c1,e0
  cx c1,f
  ry f,0.003784945883825661
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.003784945883825661
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.375, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.015139783535302645
  cx e0,f
  ry f,-0.015139783535302645
  cx c1,e0
  cx c1,f
  ry f,0.015139783535302645
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.015139783535302645
  cx c1,f
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.022413551131734478 and k=5
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 0.0009765625, controlled by [2]
  ry f,1.0944116763542225e-05
  cx c2,f
  ry f,-2.188823352708445e-05
  cx c2,f
  ry f,1.0944116763542225e-05
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 1.0, controlled by [0]
  ry f,0.011206775565867239
  cx c4,f
  ry f,-0.022413551131734478
  cx c4,f
  ry f,0.011206775565867239
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 3.0517578125e-05, controlled by [3]
  ry f,3.4200364886069455e-07
  cx c1,f
  ry f,-6.840072977213891e-07
  cx c1,f
  ry f,3.4200364886069455e-07
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 0.03125, controlled by [1]
  ry f,0.0003502117364333512
  cx c3,f
  ry f,-0.0007004234728667024
  cx c3,f
  ry f,0.0003502117364333512
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 6.5625, controlled by [0, 1]
  ry f,0.07354446465100375
  cx c3,f
  ry f,-0.07354446465100375
  cx c4,c3
  cx c4,f
  ry f,0.07354446465100375
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.07354446465100375
  cx c4,f
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 0.205078125, controlled by [1, 2]
  ry f,0.0022982645203438673
  cx c2,f
  ry f,-0.0022982645203438673
  cx c3,c2
  cx c3,f
  ry f,0.0022982645203438673
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-0.0022982645203438673
  cx c3,f
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 0.0640869140625, controlled by [1, 3]
  ry f,0.0007182076626074586
  cx c1,f
  ry f,-0.0007182076626074586
  cx c3,c1
  cx c3,f
  ry f,0.0007182076626074586
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-0.0007182076626074586
  cx c3,f
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 0.00640869140625, controlled by [2, 3]
  ry f,7.182076626074585e-05
  cx c1,f
  ry f,-7.182076626074585e-05
  cx c2,c1
  cx c2,f
  ry f,7.182076626074585e-05
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-7.182076626074585e-05
  cx c2,f
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 0.802001953125, controlled by [0, 3]
  ry f,0.008987855892059053
  cx c1,f
  ry f,-0.008987855892059053
  cx c4,c1
  cx c4,f
  ry f,0.008987855892059053
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.008987855892059053
  cx c4,f
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 2.05078125, controlled by [0, 2]
  ry f,0.022982645203438675
  cx c2,f
  ry f,-0.022982645203438675
  cx c4,c2
  cx c4,f
  ry f,0.022982645203438675
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.022982645203438675
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 6.5625, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.07354446465100375
  cx e0,f
  ry f,-0.07354446465100375
  cx c2,e0
  cx c2,f
  ry f,0.07354446465100375
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.07354446465100375
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 1.0546875, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.011819646104625603
  cx e0,f
  ry f,-0.011819646104625603
  cx c1,e0
  cx c1,f
  ry f,0.011819646104625603
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.011819646104625603
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 0.205078125, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,0.0022982645203438673
  cx e0,f
  ry f,-0.0022982645203438673
  cx c1,e0
  cx c1,f
  ry f,0.0022982645203438673
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.0022982645203438673
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 2.87109375, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.032175703284814144
  cx e0,f
  ry f,-0.032175703284814144
  cx c1,e0
  cx c1,f
  ry f,0.032175703284814144
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.032175703284814144
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.022413551131734478 and theta = 1.7578125, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,0.019699410174376007
  cx e1,f
  ry f,-0.019699410174376007
  cx c1,e1
  cx c1,f
  ry f,0.019699410174376007
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-0.019699410174376007
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.00822964594098988 and k=7
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 6.103515625e-05, controlled by [2]
  ry f,2.511488629452478e-07
  cx c2,f
  ry f,-5.022977258904956e-07
  cx c2,f
  ry f,2.511488629452478e-07
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 1.0, controlled by [0]
  ry f,0.00411482297049494
  cx c4,f
  ry f,-0.00822964594098988
  cx c4,f
  ry f,0.00411482297049494
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 4.76837158203125e-07, controlled by [3]
  ry f,1.9621004917597486e-09
  cx c1,f
  ry f,-3.924200983519497e-09
  cx c1,f
  ry f,1.9621004917597486e-09
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 0.0078125, controlled by [1]
  ry f,3.214705445699172e-05
  cx c3,f
  ry f,-6.429410891398344e-05
  cx c3,f
  ry f,3.214705445699172e-05
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 16.078125, controlled by [0, 1]
  ry f,0.06615863807248896
  cx c3,f
  ry f,-0.06615863807248896
  cx c4,c3
  cx c4,f
  ry f,0.06615863807248896
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.06615863807248896
  cx c4,f
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 0.1256103515625, controlled by [1, 2]
  ry f,0.00051686435994132
  cx c2,f
  ry f,-0.00051686435994132
  cx c3,c2
  cx c3,f
  ry f,0.00051686435994132
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-0.00051686435994132
  cx c3,f
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 0.029439926147460938, controlled by [1, 3]
  ry f,0.00012114008436124688
  cx c1,f
  ry f,-0.00012114008436124688
  cx c3,c1
  cx c3,f
  ry f,0.00012114008436124688
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-0.00012114008436124688
  cx c3,f
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 0.0009813308715820312, controlled by [2, 3]
  ry f,4.038002812041563e-06
  cx c1,f
  ry f,-4.038002812041563e-06
  cx c2,c1
  cx c2,f
  ry f,4.038002812041563e-06
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-4.038002812041563e-06
  cx c2,f
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 1.2806968688964844, controlled by [0, 3]
  ry f,0.005269840894376201
  cx c1,f
  ry f,-0.005269840894376201
  cx c4,c1
  cx c4,f
  ry f,0.005269840894376201
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.005269840894376201
  cx c4,f
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 3.768310546875, controlled by [0, 2]
  ry f,0.0155059307982396
  cx c2,f
  ry f,-0.0155059307982396
  cx c4,c2
  cx c4,f
  ry f,0.0155059307982396
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.0155059307982396
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 29.28515625, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.12050323363203347
  cx e0,f
  ry f,-0.12050323363203347
  cx c2,e0
  cx c2,f
  ry f,0.12050323363203347
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.12050323363203347
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 3.242156982421875, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.013340902025220092
  cx e0,f
  ry f,-0.013340902025220092
  cx c1,e0
  cx c1,f
  ry f,0.013340902025220092
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.013340902025220092
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 0.228790283203125, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,0.0009414315127502615
  cx e0,f
  ry f,-0.0009414315127502615
  cx c1,e0
  cx c1,f
  ry f,0.0009414315127502615
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.0009414315127502615
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 11.524749755859375, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.04742230502461611
  cx e0,f
  ry f,-0.04742230502461611
  cx c1,e0
  cx c1,f
  ry f,0.04742230502461611
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.04742230502461611
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.00822964594098988 and theta = 14.90020751953125, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,0.06131171616650862
  cx e1,f
  ry f,-0.06131171616650862
  cx c1,e1
  cx c1,f
  ry f,0.06131171616650862
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-0.06131171616650862
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.0034548126216418775 and k=9
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 3.814697265625e-06, controlled by [2]
  ry f,6.589532130512004e-09
  cx c2,f
  ry f,-1.3179064261024008e-08
  cx c2,f
  ry f,6.589532130512004e-09
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 1.0, controlled by [0]
  ry f,0.0017274063108209388
  cx c4,f
  ry f,-0.0034548126216418775
  cx c4,f
  ry f,0.0017274063108209388
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 7.450580596923828e-09, controlled by [3]
  ry f,1.2870179942406258e-11
  cx c1,f
  ry f,-2.5740359884812516e-11
  cx c1,f
  ry f,1.2870179942406258e-11
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 0.001953125, controlled by [1]
  ry f,3.373840450822146e-06
  cx c3,f
  ry f,-6.747680901644292e-06
  cx c3,f
  ry f,3.373840450822146e-06
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 37.44140625, controlled by [0, 1]
  ry f,0.06467652144226053
  cx c3,f
  ry f,-0.06467652144226053
  cx c4,c3
  cx c4,f
  ry f,0.06467652144226053
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.06467652144226053
  cx c4,f
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 0.07312774658203125, controlled by [1, 2]
  ry f,0.0001263213309419151
  cx c2,f
  ry f,-0.0001263213309419151
  cx c3,c2
  cx c3,f
  ry f,0.0001263213309419151
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-0.0001263213309419151
  cx c3,f
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 0.012598782777786255, controlled by [1, 3]
  ry f,2.1763216879010134e-05
  cx c1,f
  ry f,-2.1763216879010134e-05
  cx c3,c1
  cx c3,f
  ry f,2.1763216879010134e-05
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-2.1763216879010134e-05
  cx c3,f
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 0.00014282763004302979, controlled by [2, 3]
  ry f,2.4672134949592794e-07
  cx c1,f
  ry f,-2.4672134949592794e-07
  cx c2,c1
  cx c2,f
  ry f,2.4672134949592794e-07
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-2.4672134949592794e-07
  cx c2,f
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 1.8865075707435608, controlled by [0, 3]
  ry f,0.0032587650831139054
  cx c1,f
  ry f,-0.0032587650831139054
  cx c4,c1
  cx c4,f
  ry f,0.0032587650831139054
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.0032587650831139054
  cx c4,f
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 6.4505767822265625, controlled by [0, 2]
  ry f,0.011142767042053189
  cx c2,f
  ry f,-0.011142767042053189
  cx c4,c2
  cx c4,f
  ry f,0.011142767042053189
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.011142767042053189
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 108.9697265625, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.18823499335249458
  cx e0,f
  ry f,-0.18823499335249458
  cx c2,e0
  cx c2,f
  ry f,0.18823499335249458
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.18823499335249458
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 8.23084831237793, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.01421801931821151
  cx e0,f
  ry f,-0.01421801931821151
  cx c1,e0
  cx c1,f
  ry f,0.01421801931821151
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.01421801931821151
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 0.2128314971923828, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,0.000367646471391591
  cx e0,f
  ry f,-0.000367646471391591
  cx c1,e0
  cx c1,f
  ry f,0.000367646471391591
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.000367646471391591
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 38.66721153259277, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.06679398522324895
  cx e0,f
  ry f,-0.06679398522324895
  cx c1,e0
  cx c1,f
  ry f,0.06679398522324895
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.06679398522324895
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.0034548126216418775 and theta = 83.47841262817383, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,0.14420113679122182
  cx e1,f
  ry f,-0.14420113679122182
  cx c1,e1
  cx c1,f
  ry f,0.14420113679122182
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-0.14420113679122182
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.0015692661035767927 and k=11
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 2.384185791015625e-07, controlled by [2]
  ry f,1.8707109732351216e-10
  cx c2,f
  ry f,-3.741421946470243e-10
  cx c2,f
  ry f,1.8707109732351216e-10
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 1.0, controlled by [0]
  ry f,0.0007846330517883963
  cx c4,f
  ry f,-0.0015692661035767927
  cx c4,f
  ry f,0.0007846330517883963
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 1.1641532182693481e-10, controlled by [3]
  ry f,9.134330923999617e-14
  cx c1,f
  ry f,-1.8268661847999234e-13
  cx c1,f
  ry f,9.134330923999617e-14
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 0.00048828125, controlled by [1]
  ry f,3.831216073185529e-07
  cx c3,f
  ry f,-7.662432146371058e-07
  cx c3,f
  ry f,3.831216073185529e-07
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 85.4970703125, controlled by [0, 1]
  ry f,0.06708382719826397
  cx c3,f
  ry f,-0.06708382719826397
  cx c4,c3
  cx c4,f
  ry f,0.06708382719826397
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.06708382719826397
  cx c4,f
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 0.04174661636352539, controlled by [1, 2]
  ry f,3.275577499915233e-05
  cx c2,f
  ry f,-3.275577499915233e-05
  cx c3,c2
  cx c3,f
  ry f,3.275577499915233e-05
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-3.275577499915233e-05
  cx c3,f
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 0.00519606051966548, controlled by [1, 3]
  ry f,4.077000822822326e-06
  cx c1,f
  ry f,-4.077000822822326e-06
  cx c3,c1
  cx c3,f
  ry f,4.077000822822326e-06
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-4.077000822822326e-06
  cx c3,f
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 2.0384090021252632e-05, controlled by [2, 3]
  ry f,1.599403076130485e-08
  cx c1,f
  ry f,-1.599403076130485e-08
  cx c2,c1
  cx c2,f
  ry f,1.599403076130485e-08
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-1.599403076130485e-08
  cx c2,f
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 2.653236153535545, controlled by [0, 3]
  ry f,0.0020818167802639006
  cx c1,f
  ry f,-0.0020818167802639006
  cx c4,c1
  cx c4,f
  ry f,0.0020818167802639006
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.0020818167802639006
  cx c4,f
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 10.641531944274902, controlled by [0, 2]
  ry f,0.008349697685140123
  cx c2,f
  ry f,-0.008349697685140123
  cx c4,c2
  cx c4,f
  ry f,0.008349697685140123
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.008349697685140123
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 374.25059509277344, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.29364938656126627
  cx e0,f
  ry f,-0.29364938656126627
  cx c2,e0
  cx c2,f
  ry f,0.29364938656126627
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.29364938656126627
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 18.91986123472452, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.014845148460014878
  cx e0,f
  ry f,-0.014845148460014878
  cx c1,e0
  cx c1,f
  ry f,0.014845148460014878
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.014845148460014878
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 0.18273954838514328, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,0.0001433834895318683
  cx e0,f
  ry f,-0.0001433834895318683
  cx c1,e0
  cx c1,f
  ry f,0.0001433834895318683
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.0001433834895318683
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 119.47893822938204, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.09374712392735733
  cx e0,f
  ry f,-0.09374712392735733
  cx c1,e0
  cx c1,f
  ry f,0.09374712392735733
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.09374712392735733
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.0015692661035767927 and theta = 394.29268799722195, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,0.30937507508111023
  cx e1,f
  ry f,-0.30937507508111023
  cx c1,e1
  cx c1,f
  ry f,0.30937507508111023
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-0.30937507508111023
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.0007508222406744251 and k=13
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 1.4901161193847656e-08, controlled by [2]
  ry f,5.594061618107744e-12
  cx c2,f
  ry f,-1.1188123236215489e-11
  cx c2,f
  ry f,5.594061618107744e-12
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 1.0, controlled by [0]
  ry f,0.00037541112033721255
  cx c4,f
  ry f,-0.0007508222406744251
  cx c4,f
  ry f,0.00037541112033721255
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 1.8189894035458565e-12, controlled by [3]
  ry f,6.82868849866668e-16
  cx c1,f
  ry f,-1.365737699733336e-15
  cx c1,f
  ry f,6.82868849866668e-16
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 0.0001220703125, controlled by [1]
  ry f,4.582655277553864e-08
  cx c3,f
  ry f,-9.165310555107728e-08
  cx c3,f
  ry f,4.582655277553864e-08
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 193.619384765625, controlled by [0, 1]
  ry f,0.0726868701538651
  cx c3,f
  ry f,-0.0726868701538651
  cx c4,c3
  cx c4,f
  ry f,0.0726868701538651
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.0726868701538651
  cx c4,f
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 0.023635178804397583, controlled by [1, 2]
  ry f,8.872908954329236e-06
  cx c2,f
  ry f,-8.872908954329236e-06
  cx c3,c2
  cx c3,f
  ry f,8.872908954329236e-06
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-8.872908954329236e-06
  cx c3,f
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 0.0020983757349313237, controlled by [1, 3]
  ry f,7.8775358553899e-07
  cx c1,f
  ry f,-7.8775358553899e-07
  cx c3,c1
  cx c3,f
  ry f,7.8775358553899e-07
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-7.8775358553899e-07
  cx c3,f
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 2.8851536626461893e-06, controlled by [2, 3]
  ry f,1.083118768839018e-09
  cx c1,f
  ry f,-1.083118768839018e-09
  cx c2,c1
  cx c2,f
  ry f,1.083118768839018e-09
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-1.083118768839018e-09
  cx c2,f
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 3.623627006963943, controlled by [0, 3]
  ry f,0.0013603498743685143
  cx c1,f
  ry f,-0.0013603498743685143
  cx c4,c1
  cx c4,f
  ry f,0.0013603498743685143
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.0013603498743685143
  cx c4,f
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 17.189894020557404, controlled by [0, 2]
  ry f,0.006453277372735406
  cx c2,f
  ry f,-0.006453277372735406
  cx c4,c2
  cx c4,f
  ry f,0.006453277372735406
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.006453277372735406
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 1231.9257259368896, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.4624786169462016
  cx e0,f
  ry f,-0.4624786169462016
  cx c2,e0
  cx c2,f
  ry f,0.4624786169462016
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.4624786169462016
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 40.982923643896356, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.015385445279849568
  cx e0,f
  ry f,-0.015385445279849568
  cx c1,e0
  cx c1,f
  ry f,0.015385445279849568
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.015385445279849568
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 0.15038155834190547, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,5.645490929519063e-05
  cx e0,f
  ry f,-5.645490929519063e-05
  cx c1,e0
  cx c1,f
  ry f,5.645490929519063e-05
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-5.645490929519063e-05
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 352.6813772704918, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.13240051096318645
  cx e0,f
  ry f,-0.13240051096318645
  cx c1,e0
  cx c1,f
  ry f,0.13240051096318645
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.13240051096318645
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.0007508222406744251 and theta = 1698.9090338512324, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,0.6377893437491026
  cx e1,f
  ry f,-0.6377893437491026
  cx c1,e1
  cx c1,f
  ry f,0.6377893437491026
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-0.6377893437491026
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.00037272137648164806 and k=15
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 9.313225746154785e-10, controlled by [2]
  ry f,1.7356191597955677e-13
  cx c2,f
  ry f,-3.4712383195911353e-13
  cx c2,f
  ry f,1.7356191597955677e-13
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 1.0, controlled by [0]
  ry f,0.00018636068824082403
  cx c4,f
  ry f,-0.00037272137648164806
  cx c4,f
  ry f,0.00018636068824082403
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 2.842170943040401e-14, controlled by [3]
  ry f,5.2966893304308095e-18
  cx c1,f
  ry f,-1.0593378660861619e-17
  cx c1,f
  ry f,5.2966893304308095e-18
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 3.0517578125e-05, controlled by [1]
  ry f,5.687276862818116e-09
  cx c3,f
  ry f,-1.1374553725636232e-08
  cx c3,f
  ry f,5.687276862818116e-09
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 264.68316650390625, controlled by [0, 1]
  ry f,0.04932653707542859
  cx c3,f
  ry f,-0.04932653707542859
  cx c4,c3
  cx c4,f
  ry f,0.04932653707542859
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.04932653707542859
  cx c4,f
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 0.008077489212155342, controlled by [1, 2]
  ry f,1.505326448835101e-06
  cx c2,f
  ry f,-1.505326448835101e-06
  cx c3,c2
  cx c3,f
  ry f,1.505326448835101e-06
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-1.505326448835101e-06
  cx c3,f
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 0.0001811905661952551, controlled by [1, 3]
  ry f,3.376679861889232e-08
  cx c1,f
  ry f,-3.376679861889232e-08
  cx c3,c1
  cx c3,f
  ry f,3.376679861889232e-08
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-3.376679861889232e-08
  cx c3,f
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 2.4650540808579535e-07, controlled by [2, 3]
  ry f,4.593891750595401e-11
  cx c1,f
  ry f,-4.593891750595401e-11
  cx c2,c1
  cx c2,f
  ry f,4.593891750595401e-11
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-4.593891750595401e-11
  cx c2,f
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 2.166231035682131, controlled by [0, 3]
  ry f,0.0004037003066983549
  cx c1,f
  ry f,-0.0004037003066983549
  cx c4,c1
  cx c4,f
  ry f,0.0004037003066983549
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.0004037003066983549
  cx c4,f
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 5.937252473086119, controlled by [0, 2]
  ry f,0.0011064704571438636
  cx c2,f
  ry f,-0.0011064704571438636
  cx c4,c2
  cx c4,f
  ry f,0.0011064704571438636
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.0011064704571438636
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 3956.1822758316994, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.737276851730145
  cx e0,f
  ry f,-0.737276851730145
  cx c2,e0
  cx c2,f
  ry f,0.737276851730145
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.737276851730145
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 85.45104592183452, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.015924715728891342
  cx e0,f
  ry f,-0.015924715728891342
  cx c1,e0
  cx c1,f
  ry f,0.015924715728891342
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.015924715728891342
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 0.12073310167943418, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,2.2499903922428742e-05
  cx e0,f
  ry f,-2.2499903922428742e-05
  cx c1,e0
  cx c1,f
  ry f,2.2499903922428742e-05
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-2.2499903922428742e-05
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 1012.0440730265236, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.1886052299792697
  cx e0,f
  ry f,-0.1886052299792697
  cx c1,e0
  cx c1,f
  ry f,0.1886052299792697
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.1886052299792697
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.00037272137648164806 and theta = 6921.713237936729, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,1.2899352428275113
  cx e1,f
  ry f,-1.2899352428275113
  cx c1,e1
  cx c1,f
  ry f,1.2899352428275113
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-1.2899352428275113
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.00019018562063034358 and k=17
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 5.820766091346741e-11, controlled by [2]
  ry f,5.535130058134195e-15
  cx c2,f
  ry f,-1.107026011626839e-14
  cx c2,f
  ry f,5.535130058134195e-15
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 1.0, controlled by [0]
  ry f,9.509281031517179e-05
  cx c4,f
  ry f,-0.00019018562063034358
  cx c4,f
  ry f,9.509281031517179e-05
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 4.440892098500626e-16, controlled by [3]
  ry f,4.2229690995286524e-20
  cx c1,f
  ry f,-8.445938199057305e-20
  cx c1,f
  ry f,4.2229690995286524e-20
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 7.62939453125e-06, controlled by [1]
  ry f,7.255005669797652e-10
  cx c3,f
  ry f,-1.4510011339595305e-09
  cx c3,f
  ry f,7.255005669797652e-10
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 6768.765151977539, controlled by [0, 1]
  ry f,0.6436609006649451
  cx c3,f
  ry f,-0.6436609006649451
  cx c4,c3
  cx c4,f
  ry f,0.6436609006649451
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.6436609006649451
  cx c4,f
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 0.05164157983381301, controlled by [1, 2]
  ry f,4.9107429555125815e-06
  cx c2,f
  ry f,-4.9107429555125815e-06
  cx c3,c2
  cx c3,f
  ry f,4.9107429555125815e-06
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-4.9107429555125815e-06
  cx c3,f
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 0.00154040593483451, controlled by [1, 3]
  ry f,1.4648152936958293e-07
  cx c1,f
  ry f,-1.4648152936958293e-07
  cx c3,c1
  cx c3,f
  ry f,1.4648152936958293e-07
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-1.4648152936958293e-07
  cx c3,f
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 3.939939867692033e-07, controlled by [2, 3]
  ry f,3.746599544916215e-11
  cx c1,f
  ry f,-3.746599544916215e-11
  cx c2,c1
  cx c2,f
  ry f,3.746599544916215e-11
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-3.746599544916215e-11
  cx c2,f
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 11.897852823522332, controlled by [0, 3]
  ry f,0.00113140026170504
  cx c1,f
  ry f,-0.00113140026170504
  cx c4,c1
  cx c4,f
  ry f,0.00113140026170504
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.00113140026170504
  cx c4,f
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 201.9040866906289, controlled by [0, 2]
  ry f,0.019199627017529974
  cx c2,f
  ry f,-0.019199627017529974
  cx c4,c2
  cx c4,f
  ry f,0.019199627017529974
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.019199627017529974
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 13537.134441059083, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,1.2872841576146103
  cx e0,f
  ry f,-1.2872841576146103
  cx c2,e0
  cx c2,f
  ry f,1.2872841576146103
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-1.2872841576146103
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 155.54946850839605, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.014791636103494692
  cx e0,f
  ry f,-0.014791636103494692
  cx c1,e0
  cx c1,f
  ry f,0.014791636103494692
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.014791636103494692
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 0.10328013947341219, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,9.821198712269671e-06
  cx e0,f
  ry f,-9.821198712269671e-06
  cx c1,e0
  cx c1,f
  ry f,9.821198712269671e-06
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-9.821198712269671e-06
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 2854.6630634657704, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.2714579332078777
  cx e0,f
  ry f,-0.2714579332078777
  cx c1,e0
  cx c1,f
  ry f,0.2714579332078777
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.2714579332078777
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.00019018562063034358 and theta = 25262.5283152176, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,2.4022848131606436
  cx e1,f
  ry f,-2.4022848131606436
  cx c1,e1
  cx c1,f
  ry f,2.4022848131606436
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-2.4022848131606436
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=9.91354945717793e-05 and k=19
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 3.637978807091713e-12, controlled by [2]
  ry f,1.8032641414134432e-16
  cx c2,f
  ry f,-3.6065282828268865e-16
  cx c2,f
  ry f,1.8032641414134432e-16
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 1.0, controlled by [0]
  ry f,4.956774728588965e-05
  cx c4,f
  ry f,-9.91354945717793e-05
  cx c4,f
  ry f,4.956774728588965e-05
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 6.938893903907228e-18, controlled by [3]
  ry f,3.4394533947247376e-22
  cx c1,f
  ry f,-6.878906789449475e-22
  cx c1,f
  ry f,3.4394533947247376e-22
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 1.9073486328125e-06, controlled by [1]
  ry f,9.454297501733713e-11
  cx c3,f
  ry f,-1.8908595003467427e-10
  cx c3,f
  ry f,9.454297501733713e-11
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 1810092.6205329895, controlled by [0, 1]
  ry f,89.72221357863297
  cx c3,f
  ry f,-89.72221357863297
  cx c4,c3
  cx c4,f
  ry f,89.72221357863297
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-89.72221357863297
  cx c4,f
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 3.452477685037593, controlled by [1, 2]
  ry f,0.00017113154140211673
  cx c2,f
  ry f,-0.00017113154140211673
  cx c3,c2
  cx c3,f
  ry f,0.00017113154140211673
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-0.00017113154140211673
  cx c3,f
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 0.03215989769413105, controlled by [1, 3]
  ry f,1.5940936816427533e-06
  cx c1,f
  ry f,-1.5940936816427533e-06
  cx c3,c1
  cx c3,f
  ry f,1.5940936816427533e-06
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-1.5940936816427533e-06
  cx c3,f
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 6.585078592372118e-06, controlled by [2, 3]
  ry f,3.264075115244231e-10
  cx c1,f
  ry f,-3.264075115244231e-10
  cx c2,c1
  cx c2,f
  ry f,3.264075115244231e-10
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-3.264075115244231e-10
  cx c2,f
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 147.85109547562024, controlled by [0, 3]
  ry f,0.007328645736477487
  cx c1,f
  ry f,-0.007328645736477487
  cx c4,c1
  cx c4,f
  ry f,0.007328645736477487
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.007328645736477487
  cx c4,f
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 16861.04844226058, controlled by [0, 2]
  ry f,0.8357641881611159
  cx c2,f
  ry f,-0.8357641881611159
  cx c4,c2
  cx c4,f
  ry f,0.8357641881611159
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.8357641881611159
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 885537.0450127579, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,43.89407645948587
  cx e0,f
  ry f,-43.89407645948587
  cx c2,e0
  cx c2,f
  ry f,43.89407645948587
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-43.89407645948587
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 6059.137521211883, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.30033779742188244
  cx e0,f
  ry f,-0.30033779742188244
  cx c1,e0
  cx c1,f
  ry f,0.30033779742188244
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.30033779742188244
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 1.689027872109905, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,8.372130672356772e-05
  cx e0,f
  ry f,-8.372130672356772e-05
  cx c1,e0
  cx c1,f
  ry f,8.372130672356772e-05
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-8.372130672356772e-05
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 248847.47238567768, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,12.334808623945674
  cx e0,f
  ry f,-12.334808623945674
  cx c1,e0
  cx c1,f
  ry f,12.334808623945674
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-12.334808623945674
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 9.91354945717793e-05 and theta = 829227.6573584364, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,41.10294696241327
  cx e1,f
  ry f,-41.10294696241327
  cx c1,e1
  cx c1,f
  ry f,41.10294696241327
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-41.10294696241327
  cx c1,f
  toffoli c2,e0,e1
  toffoli c4,c3,e0
  

.thapliyal_division
  ##### un-Thapliyal division with N = [d0, ..., d4], O = [c0, ..., c4] and D = [b0, ..., b4]
  
  ### Un-iteration 5, where Y = [d0, ..., d4] and R = c0
  
  # x on Z
  x c0
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x d0
  x d1
  x d2
  x d3
  x d4
  cx b1,d1
  cx b2,d2
  cx b3,d3
  cx b4,d4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d0,b1
  toffoli b1,d1,b2
  toffoli b2,d2,b3
  toffoli b3,d3,b4
  toffoli c0,b4,d4
  toffoli b3,d3,b4
  toffoli c0,b3,d3
  toffoli b2,d2,b3
  toffoli c0,b2,d2
  toffoli b1,d1,b2
  toffoli c0,b1,d1
  toffoli b0,d0,b1
  toffoli c0,b0,d0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d4
  cx b3,d3
  cx b2,d2
  cx b1,d1
  x d0
  x d1
  x d2
  x d3
  x d4
  
  # un-SUB: Y - D, on Y
  cx b1,d1
  cx b2,d2
  cx b3,d3
  cx b4,d4
  cx b4,c0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d0,b1
  toffoli b1,d1,b2
  toffoli b2,d2,b3
  toffoli b3,d3,b4
  toffoli b4,d4,c0
  cx b4,d4
  toffoli b3,d3,b4
  cx b3,d3
  toffoli b2,d2,b3
  cx b2,d2
  toffoli b1,d1,b2
  cx b1,d1
  toffoli b0,d0,b1
  cx b0,d0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d4
  cx b3,d3
  cx b2,d2
  cx b1,d1
  
  ### Un-iteration 4, where Y = [d1, ..., c0] and R = c1
  
  # x on Z
  x c1
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x d1
  x d2
  x d3
  x d4
  x c0
  cx b1,d2
  cx b2,d3
  cx b3,d4
  cx b4,c0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d1,b1
  toffoli b1,d2,b2
  toffoli b2,d3,b3
  toffoli b3,d4,b4
  toffoli c1,b4,c0
  toffoli b3,d4,b4
  toffoli c1,b3,d4
  toffoli b2,d3,b3
  toffoli c1,b2,d3
  toffoli b1,d2,b2
  toffoli c1,b1,d2
  toffoli b0,d1,b1
  toffoli c1,b0,d1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c0
  cx b3,d4
  cx b2,d3
  cx b1,d2
  x d1
  x d2
  x d3
  x d4
  x c0
  
  # un-SUB: Y - D, on Y
  cx b1,d2
  cx b2,d3
  cx b3,d4
  cx b4,c0
  cx b4,c1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d1,b1
  toffoli b1,d2,b2
  toffoli b2,d3,b3
  toffoli b3,d4,b4
  toffoli b4,c0,c1
  cx b4,c0
  toffoli b3,d4,b4
  cx b3,d4
  toffoli b2,d3,b3
  cx b2,d3
  toffoli b1,d2,b2
  cx b1,d2
  toffoli b0,d1,b1
  cx b0,d1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c0
  cx b3,d4
  cx b2,d3
  cx b1,d2
  
  ### Un-iteration 3, where Y = [d2, ..., c1] and R = c2
  
  # x on Z
  x c2
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x d2
  x d3
  x d4
  x c0
  x c1
  cx b1,d3
  cx b2,d4
  cx b3,c0
  cx b4,c1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d2,b1
  toffoli b1,d3,b2
  toffoli b2,d4,b3
  toffoli b3,c0,b4
  toffoli c2,b4,c1
  toffoli b3,c0,b4
  toffoli c2,b3,c0
  toffoli b2,d4,b3
  toffoli c2,b2,d4
  toffoli b1,d3,b2
  toffoli c2,b1,d3
  toffoli b0,d2,b1
  toffoli c2,b0,d2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c1
  cx b3,c0
  cx b2,d4
  cx b1,d3
  x d2
  x d3
  x d4
  x c0
  x c1
  
  # un-SUB: Y - D, on Y
  cx b1,d3
  cx b2,d4
  cx b3,c0
  cx b4,c1
  cx b4,c2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d2,b1
  toffoli b1,d3,b2
  toffoli b2,d4,b3
  toffoli b3,c0,b4
  toffoli b4,c1,c2
  cx b4,c1
  toffoli b3,c0,b4
  cx b3,c0
  toffoli b2,d4,b3
  cx b2,d4
  toffoli b1,d3,b2
  cx b1,d3
  toffoli b0,d2,b1
  cx b0,d2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c1
  cx b3,c0
  cx b2,d4
  cx b1,d3
  
  ### Un-iteration 2, where Y = [d3, ..., c2] and R = c3
  
  # x on Z
  x c3
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x d3
  x d4
  x c0
  x c1
  x c2
  cx b1,d4
  cx b2,c0
  cx b3,c1
  cx b4,c2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d3,b1
  toffoli b1,d4,b2
  toffoli b2,c0,b3
  toffoli b3,c1,b4
  toffoli c3,b4,c2
  toffoli b3,c1,b4
  toffoli c3,b3,c1
  toffoli b2,c0,b3
  toffoli c3,b2,c0
  toffoli b1,d4,b2
  toffoli c3,b1,d4
  toffoli b0,d3,b1
  toffoli c3,b0,d3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c2
  cx b3,c1
  cx b2,c0
  cx b1,d4
  x d3
  x d4
  x c0
  x c1
  x c2
  
  # un-SUB: Y - D, on Y
  cx b1,d4
  cx b2,c0
  cx b3,c1
  cx b4,c2
  cx b4,c3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d3,b1
  toffoli b1,d4,b2
  toffoli b2,c0,b3
  toffoli b3,c1,b4
  toffoli b4,c2,c3
  cx b4,c2
  toffoli b3,c1,b4
  cx b3,c1
  toffoli b2,c0,b3
  cx b2,c0
  toffoli b1,d4,b2
  cx b1,d4
  toffoli b0,d3,b1
  cx b0,d3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c2
  cx b3,c1
  cx b2,c0
  cx b1,d4
  
  ### Un-iteration 1, where Y = [d4, ..., c3] and R = c4
  
  # x on Z
  x c4
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x d4
  x c0
  x c1
  x c2
  x c3
  cx b1,c0
  cx b2,c1
  cx b3,c2
  cx b4,c3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d4,b1
  toffoli b1,c0,b2
  toffoli b2,c1,b3
  toffoli b3,c2,b4
  toffoli c4,b4,c3
  toffoli b3,c2,b4
  toffoli c4,b3,c2
  toffoli b2,c1,b3
  toffoli c4,b2,c1
  toffoli b1,c0,b2
  toffoli c4,b1,c0
  toffoli b0,d4,b1
  toffoli c4,b0,d4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c3
  cx b3,c2
  cx b2,c1
  cx b1,c0
  x d4
  x c0
  x c1
  x c2
  x c3
  
  # un-SUB: Y - D, on Y
  cx b1,c0
  cx b2,c1
  cx b3,c2
  cx b4,c3
  cx b4,c4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,d4,b1
  toffoli b1,c0,b2
  toffoli b2,c1,b3
  toffoli b3,c2,b4
  toffoli b4,c3,c4
  cx b4,c3
  toffoli b3,c2,b4
  cx b3,c2
  toffoli b2,c1,b3
  cx b2,c1
  toffoli b1,c0,b2
  cx b1,c0
  toffoli b0,d4,b1
  cx b0,d4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c3
  cx b3,c2
  cx b2,c1
  cx b1,c0

.reverse
  swap b0,b3
  swap b1,b2

.qft
  h b0
  cr b1,b0
  cr b2,b0
  cr b3,b0
  h b1
  cr b2,b1
  cr b3,b1
  h b2
  cr b3,b2
  h b3

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/2) on a0 and a1, controlled by b3
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b3,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b3,a0
  tdag b3
  tdag a1
  cx b3,a1
  t a1
  cx b3,a1
  
  # Toffoli(a,b,c)
  toffoli b3,a0,a1
  
  # cRx(a,b,4.712)
  rx a0,1.178
  cz b3,a0
  rx a0,-2.356
  cz b3,a0
  rx a0,1.178
  
  # Toffoli(a,b,c)
  toffoli b3,a0,a1
  
  # Rz(a,-15.064)
  rz b3,-15.064
  
  # cRx(a,b,7.856)
  rx a0,1.964
  cz b3,a0
  rx a0,-3.928
  cz b3,a0
  rx a0,1.964
  
  # Rz(a,-3.0)
  rz b3,-3.0
  
  # cRx(a,c,1.568)
  rx a1,0.392
  cz b3,a1
  rx a1,-0.784
  cz b3,a1
  rx a1,0.392
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b3,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b3,a0
  tdag b3
  tdag a1
  cx b3,a1
  t a1
  cx b3,a1

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/4) on a0 and a1, controlled by b2
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b2,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b2,a0
  tdag b2
  tdag a1
  cx b2,a1
  t a1
  cx b2,a1
  
  # Toffoli(a,b,c)
  toffoli b2,a0,a1
  
  # cRx(a,b,2.356)
  rx a0,0.589
  cz b2,a0
  rx a0,-1.178
  cz b2,a0
  rx a0,0.589
  
  # Toffoli(a,b,c)
  toffoli b2,a0,a1
  
  # Rz(a,-7.532)
  rz b2,-7.532
  
  # cRx(a,b,3.928)
  rx a0,0.982
  cz b2,a0
  rx a0,-1.964
  cz b2,a0
  rx a0,0.982
  
  # Rz(a,-1.5)
  rz b2,-1.5
  
  # cRx(a,c,0.784)
  rx a1,0.196
  cz b2,a1
  rx a1,-0.392
  cz b2,a1
  rx a1,0.196
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b2,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b2,a0
  tdag b2
  tdag a1
  cx b2,a1
  t a1
  cx b2,a1

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/8) on a0 and a1, controlled by b1
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b1,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b1,a0
  tdag b1
  tdag a1
  cx b1,a1
  t a1
  cx b1,a1
  
  # Toffoli(a,b,c)
  toffoli b1,a0,a1
  
  # cRx(a,b,1.178)
  rx a0,0.2945
  cz b1,a0
  rx a0,-0.589
  cz b1,a0
  rx a0,0.2945
  
  # Toffoli(a,b,c)
  toffoli b1,a0,a1
  
  # Rz(a,-3.766)
  rz b1,-3.766
  
  # cRx(a,b,1.964)
  rx a0,0.491
  cz b1,a0
  rx a0,-0.982
  cz b1,a0
  rx a0,0.491
  
  # Rz(a,-0.75)
  rz b1,-0.75
  
  # cX(a,c)
  cx b1,a1
  
  # cRx(a,c,0.392)
  rx a1,0.098
  cz b1,a1
  rx a1,-0.196
  cz b1,a1
  rx a1,0.098
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b1,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b1,a0
  tdag b1
  tdag a1
  cx b1,a1
  t a1
  cx b1,a1

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/16) on a0 and a1, controlled by b0
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b0,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b0,a0
  tdag b0
  tdag a1
  cx b0,a1
  t a1
  cx b0,a1
  
  # Toffoli(a,b,c)
  toffoli b0,a0,a1
  
  # cRx(a,b,0.589)
  rx a0,0.14725
  cz b0,a0
  rx a0,-0.2945
  cz b0,a0
  rx a0,0.14725
  
  # Toffoli(a,b,c)
  toffoli b0,a0,a1
  
  # Rz(a,-1.883)
  rz b0,-1.883
  
  # cRx(a,b,0.982)
  rx a0,0.2455
  cz b0,a0
  rx a0,-0.491
  cz b0,a0
  rx a0,0.2455
  
  # Rz(a,-0.375)
  rz b0,-0.375
  
  # cV(a,c)
  t b0
  x b0
  tdag b0
  x b0
  h a1
  tdag a1
  tdag a1
  cx b0,a1
  h a1
  t a1
  cx b0,a1
  t a1
  h a1
  s a1
  x a1
  t b0
  t a1
  cx b0,a1
  tdag a1
  cx b0,a1
  t b0
  x b0
  tdag b0
  x b0
  h a1
  tdag a1
  tdag a1
  cx b0,a1
  h a1
  t a1
  cx b0,a1
  t a1
  h a1
  s a1
  x a1
  
  # cRx(a,c,0.196)
  rx a1,0.049
  cz b0,a1
  rx a1,-0.098
  cz b0,a1
  rx a1,0.049
  
  # ccZ(a,b,c)
  tdag a0
  cx a0,a1
  t a1
  cx a0,a1
  cx b0,a0
  t a0
  cx a0,a1
  tdag a1
  cx a0,a1
  cx b0,a0
  tdag b0
  tdag a1
  cx b0,a1
  t a1
  cx b0,a1

.uninit
  h b0
  h b1
  h b2
  h b3
  x d4

.result
  display
  measure
  display