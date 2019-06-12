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
  ####### Ancilla rotation using the first m terms of arcsin(x) and c=0.39269908169872414
  
  ##### Performs Ry(c*(x^k)) for c=0.39269908169872414 and k=1
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.39269908169872414 and theta = 0.25, controlled by [2]
  ry f,0.04908738521234052
  cx c2,f
  ry f,-0.09817477042468103
  cx c2,f
  ry f,0.04908738521234052
  
  # Rotation c*theta with c = 0.39269908169872414 and theta = 1.0, controlled by [0]
  ry f,0.19634954084936207
  cx c4,f
  ry f,-0.39269908169872414
  cx c4,f
  ry f,0.19634954084936207
  
  # Rotation c*theta with c = 0.39269908169872414 and theta = 0.125, controlled by [3]
  ry f,0.02454369260617026
  cx c1,f
  ry f,-0.04908738521234052
  cx c1,f
  ry f,0.02454369260617026
  
  # Rotation c*theta with c = 0.39269908169872414 and theta = 0.5, controlled by [1]
  ry f,0.09817477042468103
  cx c3,f
  ry f,-0.19634954084936207
  cx c3,f
  ry f,0.09817477042468103
  
  ##### Performs Ry(c*(x^k)) for c=0.010093189023535096 and k=3
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.015625, controlled by [2]
  ry f,7.885303924636794e-05
  cx c2,f
  ry f,-0.00015770607849273588
  cx c2,f
  ry f,7.885303924636794e-05
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 1.0, controlled by [0]
  ry f,0.005046594511767548
  cx c4,f
  ry f,-0.010093189023535096
  cx c4,f
  ry f,0.005046594511767548
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.001953125, controlled by [3]
  ry f,9.856629905795993e-06
  cx c1,f
  ry f,-1.9713259811591985e-05
  cx c1,f
  ry f,9.856629905795993e-06
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.125, controlled by [1]
  ry f,0.0006308243139709435
  cx c3,f
  ry f,-0.001261648627941887
  cx c3,f
  ry f,0.0006308243139709435
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 2.25, controlled by [0, 1]
  ry f,0.011354837651476983
  cx c3,f
  ry f,-0.011354837651476983
  cx c4,c3
  cx c4,f
  ry f,0.011354837651476983
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.011354837651476983
  cx c4,f
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.28125, controlled by [1, 2]
  ry f,0.0014193547064346229
  cx c2,f
  ry f,-0.0014193547064346229
  cx c3,c2
  cx c3,f
  ry f,0.0014193547064346229
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-0.0014193547064346229
  cx c3,f
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.1171875, controlled by [1, 3]
  ry f,0.0005913977943477596
  cx c1,f
  ry f,-0.0005913977943477596
  cx c3,c1
  cx c3,f
  ry f,0.0005913977943477596
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-0.0005913977943477596
  cx c3,f
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.03515625, controlled by [2, 3]
  ry f,0.00017741933830432786
  cx c1,f
  ry f,-0.00017741933830432786
  cx c2,c1
  cx c2,f
  ry f,0.00017741933830432786
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-0.00017741933830432786
  cx c2,f
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.421875, controlled by [0, 3]
  ry f,0.0021290320596519343
  cx c1,f
  ry f,-0.0021290320596519343
  cx c4,c1
  cx c4,f
  ry f,0.0021290320596519343
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.0021290320596519343
  cx c4,f
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.9375, controlled by [0, 2]
  ry f,0.0047311823547820766
  cx c2,f
  ry f,-0.0047311823547820766
  cx c4,c2
  cx c4,f
  ry f,0.0047311823547820766
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.0047311823547820766
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.75, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.003784945883825661
  cx e0,f
  ry f,-0.003784945883825661
  cx c2,e0
  cx c2,f
  ry f,0.003784945883825661
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.003784945883825661
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.1875, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.0009462364709564153
  cx e0,f
  ry f,-0.0009462364709564153
  cx c1,e0
  cx c1,f
  ry f,0.0009462364709564153
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.0009462364709564153
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.09375, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,0.00047311823547820764
  cx e0,f
  ry f,-0.00047311823547820764
  cx c1,e0
  cx c1,f
  ry f,0.00047311823547820764
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.00047311823547820764
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.010093189023535096 and theta = 0.375, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.0018924729419128306
  cx e0,f
  ry f,-0.0018924729419128306
  cx c1,e0
  cx c1,f
  ry f,0.0018924729419128306
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.0018924729419128306
  cx c1,f
  toffoli c4,c3,e0
  
  ##### Performs Ry(c*(x^k)) for c=0.0007004234728667024 and k=5
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 0.0009765625, controlled by [2]
  ry f,3.4200364886069455e-07
  cx c2,f
  ry f,-6.840072977213891e-07
  cx c2,f
  ry f,3.4200364886069455e-07
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 1.0, controlled by [0]
  ry f,0.0003502117364333512
  cx c4,f
  ry f,-0.0007004234728667024
  cx c4,f
  ry f,0.0003502117364333512
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 3.0517578125e-05, controlled by [3]
  ry f,1.0687614026896705e-08
  cx c1,f
  ry f,-2.137522805379341e-08
  cx c1,f
  ry f,1.0687614026896705e-08
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 0.03125, controlled by [1]
  ry f,1.0944116763542225e-05
  cx c3,f
  ry f,-2.188823352708445e-05
  cx c3,f
  ry f,1.0944116763542225e-05
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 6.5625, controlled by [0, 1]
  ry f,0.0022982645203438673
  cx c3,f
  ry f,-0.0022982645203438673
  cx c4,c3
  cx c4,f
  ry f,0.0022982645203438673
  cx c3,f
  cx c4,c3
  cx c4,f
  ry f,-0.0022982645203438673
  cx c4,f
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 0.205078125, controlled by [1, 2]
  ry f,7.182076626074585e-05
  cx c2,f
  ry f,-7.182076626074585e-05
  cx c3,c2
  cx c3,f
  ry f,7.182076626074585e-05
  cx c2,f
  cx c3,c2
  cx c3,f
  ry f,-7.182076626074585e-05
  cx c3,f
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 0.0640869140625, controlled by [1, 3]
  ry f,2.244398945648308e-05
  cx c1,f
  ry f,-2.244398945648308e-05
  cx c3,c1
  cx c3,f
  ry f,2.244398945648308e-05
  cx c1,f
  cx c3,c1
  cx c3,f
  ry f,-2.244398945648308e-05
  cx c3,f
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 0.00640869140625, controlled by [2, 3]
  ry f,2.244398945648308e-06
  cx c1,f
  ry f,-2.244398945648308e-06
  cx c2,c1
  cx c2,f
  ry f,2.244398945648308e-06
  cx c1,f
  cx c2,c1
  cx c2,f
  ry f,-2.244398945648308e-06
  cx c2,f
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 0.802001953125, controlled by [0, 3]
  ry f,0.0002808704966268454
  cx c1,f
  ry f,-0.0002808704966268454
  cx c4,c1
  cx c4,f
  ry f,0.0002808704966268454
  cx c1,f
  cx c4,c1
  cx c4,f
  ry f,-0.0002808704966268454
  cx c4,f
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 2.05078125, controlled by [0, 2]
  ry f,0.0007182076626074586
  cx c2,f
  ry f,-0.0007182076626074586
  cx c4,c2
  cx c4,f
  ry f,0.0007182076626074586
  cx c2,f
  cx c4,c2
  cx c4,f
  ry f,-0.0007182076626074586
  cx c4,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 6.5625, controlled by [0, 1, 2]
  toffoli c4,c3,e0
  ry f,0.0022982645203438673
  cx e0,f
  ry f,-0.0022982645203438673
  cx c2,e0
  cx c2,f
  ry f,0.0022982645203438673
  cx e0,f
  cx c2,e0
  cx c2,f
  ry f,-0.0022982645203438673
  cx c2,f
  toffoli c4,c3,e0
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 1.0546875, controlled by [0, 2, 3]
  toffoli c4,c2,e0
  ry f,0.0003693639407695501
  cx e0,f
  ry f,-0.0003693639407695501
  cx c1,e0
  cx c1,f
  ry f,0.0003693639407695501
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.0003693639407695501
  cx c1,f
  toffoli c4,c2,e0
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 0.205078125, controlled by [1, 2, 3]
  toffoli c3,c2,e0
  ry f,7.182076626074585e-05
  cx e0,f
  ry f,-7.182076626074585e-05
  cx c1,e0
  cx c1,f
  ry f,7.182076626074585e-05
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-7.182076626074585e-05
  cx c1,f
  toffoli c3,c2,e0
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 2.87109375, controlled by [0, 1, 3]
  toffoli c4,c3,e0
  ry f,0.001005490727650442
  cx e0,f
  ry f,-0.001005490727650442
  cx c1,e0
  cx c1,f
  ry f,0.001005490727650442
  cx e0,f
  cx c1,e0
  cx c1,f
  ry f,-0.001005490727650442
  cx c1,f
  toffoli c4,c3,e0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.0007004234728667024 and theta = 1.7578125, controlled by [0, 1, 2, 3]
  toffoli c4,c3,e0
  toffoli c2,e0,e1
  ry f,0.0006156065679492502
  cx e1,f
  ry f,-0.0006156065679492502
  cx c1,e1
  cx c1,f
  ry f,0.0006156065679492502
  cx e1,f
  cx c1,e1
  cx c1,f
  ry f,-0.0006156065679492502
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