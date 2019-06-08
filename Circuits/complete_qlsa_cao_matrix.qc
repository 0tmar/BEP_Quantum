#function: Complete QLSA for the Cao matrix

qubits 20

.init
  map q0,a0
  map q1,a1
  map q2,b0
  map q3,b1
  map q4,b2
  map q5,b3
  map q6,b4
  map q7,c0
  map q8,c1
  map q9,c2
  map q10,c3
  map q11,c4
  map q12,c5
  map q13,d0
  map q14,d1
  map q15,d2
  map q16,d3
  map q17,d4
  map q18,e
  map q19,f
  h b0
  h b1
  h b2
  h b3
  h a0
  h a1
  x c5
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
  ##### Thapliyal division with Q = [c0, ..., c5], R = [d0, ..., d4] and D = [b0, ..., b4]
  
  ### Iteration 1, where Y = [c5, ..., d3] and Z = d4
  
  # SUB: Y - D, on Y
  x c5
  x d0
  x d1
  x d2
  x d3
  cx b1,d0
  cx b2,d1
  cx b3,d2
  cx b4,d3
  cx b4,d4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c5,b1
  toffoli b1,d0,b2
  toffoli b2,d1,b3
  toffoli b3,d2,b4
  toffoli b4,d3,d4
  cx b4,d3
  toffoli b3,d2,b4
  cx b3,d2
  toffoli b2,d1,b3
  cx b2,d1
  toffoli b1,d0,b2
  cx b1,d0
  toffoli b0,c5,b1
  cx b0,c5
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d3
  cx b3,d2
  cx b2,d1
  cx b1,d0
  x c5
  x d0
  x d1
  x d2
  x d3
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,d0
  cx b2,d1
  cx b3,d2
  cx b4,d3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c5,b1
  toffoli b1,d0,b2
  toffoli b2,d1,b3
  toffoli b3,d2,b4
  toffoli d4,b4,d3
  toffoli b3,d2,b4
  toffoli d4,b3,d2
  toffoli b2,d1,b3
  toffoli d4,b2,d1
  toffoli b1,d0,b2
  toffoli d4,b1,d0
  toffoli b0,c5,b1
  toffoli d4,b0,c5
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d3
  cx b3,d2
  cx b2,d1
  cx b1,d0
  
  # x on Z
  x d4
  
  ### Iteration 2, where Y = [c4, ..., d2] and Z = d3
  
  # SUB: Y - D, on Y
  x c4
  x c5
  x d0
  x d1
  x d2
  cx b1,c5
  cx b2,d0
  cx b3,d1
  cx b4,d2
  cx b4,d3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c4,b1
  toffoli b1,c5,b2
  toffoli b2,d0,b3
  toffoli b3,d1,b4
  toffoli b4,d2,d3
  cx b4,d2
  toffoli b3,d1,b4
  cx b3,d1
  toffoli b2,d0,b3
  cx b2,d0
  toffoli b1,c5,b2
  cx b1,c5
  toffoli b0,c4,b1
  cx b0,c4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d2
  cx b3,d1
  cx b2,d0
  cx b1,c5
  x c4
  x c5
  x d0
  x d1
  x d2
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,c5
  cx b2,d0
  cx b3,d1
  cx b4,d2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c4,b1
  toffoli b1,c5,b2
  toffoli b2,d0,b3
  toffoli b3,d1,b4
  toffoli d3,b4,d2
  toffoli b3,d1,b4
  toffoli d3,b3,d1
  toffoli b2,d0,b3
  toffoli d3,b2,d0
  toffoli b1,c5,b2
  toffoli d3,b1,c5
  toffoli b0,c4,b1
  toffoli d3,b0,c4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d2
  cx b3,d1
  cx b2,d0
  cx b1,c5
  
  # x on Z
  x d3
  
  ### Iteration 3, where Y = [c3, ..., d1] and Z = d2
  
  # SUB: Y - D, on Y
  x c3
  x c4
  x c5
  x d0
  x d1
  cx b1,c4
  cx b2,c5
  cx b3,d0
  cx b4,d1
  cx b4,d2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c3,b1
  toffoli b1,c4,b2
  toffoli b2,c5,b3
  toffoli b3,d0,b4
  toffoli b4,d1,d2
  cx b4,d1
  toffoli b3,d0,b4
  cx b3,d0
  toffoli b2,c5,b3
  cx b2,c5
  toffoli b1,c4,b2
  cx b1,c4
  toffoli b0,c3,b1
  cx b0,c3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d1
  cx b3,d0
  cx b2,c5
  cx b1,c4
  x c3
  x c4
  x c5
  x d0
  x d1
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,c4
  cx b2,c5
  cx b3,d0
  cx b4,d1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c3,b1
  toffoli b1,c4,b2
  toffoli b2,c5,b3
  toffoli b3,d0,b4
  toffoli d2,b4,d1
  toffoli b3,d0,b4
  toffoli d2,b3,d0
  toffoli b2,c5,b3
  toffoli d2,b2,c5
  toffoli b1,c4,b2
  toffoli d2,b1,c4
  toffoli b0,c3,b1
  toffoli d2,b0,c3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d1
  cx b3,d0
  cx b2,c5
  cx b1,c4
  
  # x on Z
  x d2
  
  ### Iteration 4, where Y = [c2, ..., d0] and Z = d1
  
  # SUB: Y - D, on Y
  x c2
  x c3
  x c4
  x c5
  x d0
  cx b1,c3
  cx b2,c4
  cx b3,c5
  cx b4,d0
  cx b4,d1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c2,b1
  toffoli b1,c3,b2
  toffoli b2,c4,b3
  toffoli b3,c5,b4
  toffoli b4,d0,d1
  cx b4,d0
  toffoli b3,c5,b4
  cx b3,c5
  toffoli b2,c4,b3
  cx b2,c4
  toffoli b1,c3,b2
  cx b1,c3
  toffoli b0,c2,b1
  cx b0,c2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d0
  cx b3,c5
  cx b2,c4
  cx b1,c3
  x c2
  x c3
  x c4
  x c5
  x d0
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,c3
  cx b2,c4
  cx b3,c5
  cx b4,d0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c2,b1
  toffoli b1,c3,b2
  toffoli b2,c4,b3
  toffoli b3,c5,b4
  toffoli d1,b4,d0
  toffoli b3,c5,b4
  toffoli d1,b3,c5
  toffoli b2,c4,b3
  toffoli d1,b2,c4
  toffoli b1,c3,b2
  toffoli d1,b1,c3
  toffoli b0,c2,b1
  toffoli d1,b0,c2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d0
  cx b3,c5
  cx b2,c4
  cx b1,c3
  
  # x on Z
  x d1
  
  ### Iteration 5, where Y = [c1, ..., c5] and Z = d0
  
  # SUB: Y - D, on Y
  x c1
  x c2
  x c3
  x c4
  x c5
  cx b1,c2
  cx b2,c3
  cx b3,c4
  cx b4,c5
  cx b4,d0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c1,b1
  toffoli b1,c2,b2
  toffoli b2,c3,b3
  toffoli b3,c4,b4
  toffoli b4,c5,d0
  cx b4,c5
  toffoli b3,c4,b4
  cx b3,c4
  toffoli b2,c3,b3
  cx b2,c3
  toffoli b1,c2,b2
  cx b1,c2
  toffoli b0,c1,b1
  cx b0,c1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c5
  cx b3,c4
  cx b2,c3
  cx b1,c2
  x c1
  x c2
  x c3
  x c4
  x c5
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,c2
  cx b2,c3
  cx b3,c4
  cx b4,c5
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c1,b1
  toffoli b1,c2,b2
  toffoli b2,c3,b3
  toffoli b3,c4,b4
  toffoli d0,b4,c5
  toffoli b3,c4,b4
  toffoli d0,b3,c4
  toffoli b2,c3,b3
  toffoli d0,b2,c3
  toffoli b1,c2,b2
  toffoli d0,b1,c2
  toffoli b0,c1,b1
  toffoli d0,b0,c1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c5
  cx b3,c4
  cx b2,c3
  cx b1,c2
  
  # x on Z
  x d0
  
  ### Iteration 6, where Y = [c0, ..., c4] and Z = c5
  
  # SUB: Y - D, on Y
  x c0
  x c1
  x c2
  x c3
  x c4
  cx b1,c1
  cx b2,c2
  cx b3,c3
  cx b4,c4
  cx b4,c5
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c0,b1
  toffoli b1,c1,b2
  toffoli b2,c2,b3
  toffoli b3,c3,b4
  toffoli b4,c4,c5
  cx b4,c4
  toffoli b3,c3,b4
  cx b3,c3
  toffoli b2,c2,b3
  cx b2,c2
  toffoli b1,c1,b2
  cx b1,c1
  toffoli b0,c0,b1
  cx b0,c0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c4
  cx b3,c3
  cx b2,c2
  cx b1,c1
  x c0
  x c1
  x c2
  x c3
  x c4
  
  # cADD: Y (+ D), on Y, ctrl'd by Z
  cx b1,c1
  cx b2,c2
  cx b3,c3
  cx b4,c4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c0,b1
  toffoli b1,c1,b2
  toffoli b2,c2,b3
  toffoli b3,c3,b4
  toffoli c5,b4,c4
  toffoli b3,c3,b4
  toffoli c5,b3,c3
  toffoli b2,c2,b3
  toffoli c5,b2,c2
  toffoli b1,c1,b2
  toffoli c5,b1,c1
  toffoli b0,c0,b1
  toffoli c5,b0,c0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c4
  cx b3,c3
  cx b2,c2
  cx b1,c1
  
  # x on Z
  x c5

.display
  display

.ancilla_rotation
  ####### Ancilla rotation using the first m terms of arcsin(x) and c=0.7853981633974483
  
  ##### Performs Ry(c*(x^k)) for c=0.7853981633974483 and k=1
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 1.0, controlled by [0]
  ry f,0.39269908169872414
  cx d4,f
  ry f,-0.7853981633974483
  cx d4,f
  ry f,0.39269908169872414
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.5, controlled by [1]
  ry f,0.19634954084936207
  cx d3,f
  ry f,-0.39269908169872414
  cx d3,f
  ry f,0.19634954084936207
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.25, controlled by [2]
  ry f,0.09817477042468103
  cx d2,f
  ry f,-0.19634954084936207
  cx d2,f
  ry f,0.09817477042468103
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.125, controlled by [3]
  ry f,0.04908738521234052
  cx d1,f
  ry f,-0.09817477042468103
  cx d1,f
  ry f,0.04908738521234052
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.0625, controlled by [4]
  ry f,0.02454369260617026
  cx d0,f
  ry f,-0.04908738521234052
  cx d0,f
  ry f,0.02454369260617026
  
  # Rotation c*theta with c = 0.7853981633974483 and theta = 0.03125, controlled by [5]
  ry f,0.01227184630308513
  cx c5,f
  ry f,-0.02454369260617026
  cx c5,f
  ry f,0.01227184630308513
  
  ##### Performs Ry(c*(x^k)) for c=0.08074551218828077 and k=3
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 1.0, controlled by [0]
  ry f,0.040372756094140386
  cx d4,f
  ry f,-0.08074551218828077
  cx d4,f
  ry f,0.040372756094140386
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.125, controlled by [1]
  ry f,0.005046594511767548
  cx d3,f
  ry f,-0.010093189023535096
  cx d3,f
  ry f,0.005046594511767548
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.015625, controlled by [2]
  ry f,0.0006308243139709435
  cx d2,f
  ry f,-0.001261648627941887
  cx d2,f
  ry f,0.0006308243139709435
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.001953125, controlled by [3]
  ry f,7.885303924636794e-05
  cx d1,f
  ry f,-0.00015770607849273588
  cx d1,f
  ry f,7.885303924636794e-05
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.000244140625, controlled by [4]
  ry f,9.856629905795993e-06
  cx d0,f
  ry f,-1.9713259811591985e-05
  cx d0,f
  ry f,9.856629905795993e-06
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 3.0517578125e-05, controlled by [5]
  ry f,1.232078738224499e-06
  cx c5,f
  ry f,-2.464157476448998e-06
  cx c5,f
  ry f,1.232078738224499e-06
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 2.25, controlled by [0, 1]
  ry f,0.09083870121181586
  cx d3,f
  ry f,-0.09083870121181586
  cx d4,d3
  cx d4,f
  ry f,0.09083870121181586
  cx d3,f
  cx d4,d3
  cx d4,f
  ry f,-0.09083870121181586
  cx d4,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.28125, controlled by [1, 2]
  ry f,0.011354837651476983
  cx d2,f
  ry f,-0.011354837651476983
  cx d3,d2
  cx d3,f
  ry f,0.011354837651476983
  cx d2,f
  cx d3,d2
  cx d3,f
  ry f,-0.011354837651476983
  cx d3,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.1171875, controlled by [1, 3]
  ry f,0.0047311823547820766
  cx d1,f
  ry f,-0.0047311823547820766
  cx d3,d1
  cx d3,f
  ry f,0.0047311823547820766
  cx d1,f
  cx d3,d1
  cx d3,f
  ry f,-0.0047311823547820766
  cx d3,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.00054931640625, controlled by [4, 5]
  ry f,2.2177417288040983e-05
  cx c5,f
  ry f,-2.2177417288040983e-05
  cx d0,c5
  cx d0,f
  ry f,2.2177417288040983e-05
  cx c5,f
  cx d0,c5
  cx d0,f
  ry f,-2.2177417288040983e-05
  cx d0,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.052734375, controlled by [1, 4]
  ry f,0.0021290320596519343
  cx d0,f
  ry f,-0.0021290320596519343
  cx d3,d0
  cx d3,f
  ry f,0.0021290320596519343
  cx d0,f
  cx d3,d0
  cx d3,f
  ry f,-0.0021290320596519343
  cx d3,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.02490234375, controlled by [1, 5]
  ry f,0.0010053762503911912
  cx c5,f
  ry f,-0.0010053762503911912
  cx d3,c5
  cx d3,f
  ry f,0.0010053762503911912
  cx c5,f
  cx d3,c5
  cx d3,f
  ry f,-0.0010053762503911912
  cx d3,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.0146484375, controlled by [2, 4]
  ry f,0.0005913977943477596
  cx d0,f
  ry f,-0.0005913977943477596
  cx d2,d0
  cx d2,f
  ry f,0.0005913977943477596
  cx d0,f
  cx d2,d0
  cx d2,f
  ry f,-0.0005913977943477596
  cx d2,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.0966796875, controlled by [0, 5]
  ry f,0.003903225442695213
  cx c5,f
  ry f,-0.003903225442695213
  cx d4,c5
  cx d4,f
  ry f,0.003903225442695213
  cx c5,f
  cx d4,c5
  cx d4,f
  ry f,-0.003903225442695213
  cx d4,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.03515625, controlled by [2, 3]
  ry f,0.0014193547064346229
  cx d1,f
  ry f,-0.0014193547064346229
  cx d2,d1
  cx d2,f
  ry f,0.0014193547064346229
  cx d1,f
  cx d2,d1
  cx d2,f
  ry f,-0.0014193547064346229
  cx d2,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.006591796875, controlled by [2, 5]
  ry f,0.0002661290074564918
  cx c5,f
  ry f,-0.0002661290074564918
  cx d2,c5
  cx d2,f
  ry f,0.0002661290074564918
  cx c5,f
  cx d2,c5
  cx d2,f
  ry f,-0.0002661290074564918
  cx d2,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.19921875, controlled by [0, 4]
  ry f,0.00804301000312953
  cx d0,f
  ry f,-0.00804301000312953
  cx d4,d0
  cx d4,f
  ry f,0.00804301000312953
  cx d0,f
  cx d4,d0
  cx d4,f
  ry f,-0.00804301000312953
  cx d4,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.421875, controlled by [0, 3]
  ry f,0.017032256477215475
  cx d1,f
  ry f,-0.017032256477215475
  cx d4,d1
  cx d4,f
  ry f,0.017032256477215475
  cx d1,f
  cx d4,d1
  cx d4,f
  ry f,-0.017032256477215475
  cx d4,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.00439453125, controlled by [3, 4]
  ry f,0.00017741933830432786
  cx d0,f
  ry f,-0.00017741933830432786
  cx d1,d0
  cx d1,f
  ry f,0.00017741933830432786
  cx d0,f
  cx d1,d0
  cx d1,f
  ry f,-0.00017741933830432786
  cx d1,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.9375, controlled by [0, 2]
  ry f,0.03784945883825661
  cx d2,f
  ry f,-0.03784945883825661
  cx d4,d2
  cx d4,f
  ry f,0.03784945883825661
  cx d2,f
  cx d4,d2
  cx d4,f
  ry f,-0.03784945883825661
  cx d4,f
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.0018310546875, controlled by [3, 5]
  ry f,7.392472429346995e-05
  cx c5,f
  ry f,-7.392472429346995e-05
  cx d1,c5
  cx d1,f
  ry f,7.392472429346995e-05
  cx c5,f
  cx d1,c5
  cx d1,f
  ry f,-7.392472429346995e-05
  cx d1,f
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.005859375, controlled by [2, 3, 5]
  toffoli d2,d1,e
  ry f,0.00023655911773910382
  cx e,f
  ry f,-0.00023655911773910382
  cx c5,e
  cx c5,f
  ry f,0.00023655911773910382
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.00023655911773910382
  cx c5,f
  toffoli d2,d1,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.00146484375, controlled by [3, 4, 5]
  toffoli d1,d0,e
  ry f,5.9139779434775956e-05
  cx e,f
  ry f,-5.9139779434775956e-05
  cx c5,e
  cx c5,f
  ry f,5.9139779434775956e-05
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-5.9139779434775956e-05
  cx c5,f
  toffoli d1,d0,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.046875, controlled by [0, 3, 4]
  toffoli d4,d1,e
  ry f,0.0018924729419128306
  cx e,f
  ry f,-0.0018924729419128306
  cx d0,e
  cx d0,f
  ry f,0.0018924729419128306
  cx e,f
  cx d0,e
  cx d0,f
  ry f,-0.0018924729419128306
  cx d0,f
  toffoli d4,d1,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.01171875, controlled by [0, 4, 5]
  toffoli d4,d0,e
  ry f,0.00047311823547820764
  cx e,f
  ry f,-0.00047311823547820764
  cx c5,e
  cx c5,f
  ry f,0.00047311823547820764
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.00047311823547820764
  cx c5,f
  toffoli d4,d0,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.75, controlled by [0, 1, 2]
  toffoli d4,d3,e
  ry f,0.03027956707060529
  cx e,f
  ry f,-0.03027956707060529
  cx d2,e
  cx d2,f
  ry f,0.03027956707060529
  cx e,f
  cx d2,e
  cx d2,f
  ry f,-0.03027956707060529
  cx d2,f
  toffoli d4,d3,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.1875, controlled by [0, 2, 3]
  toffoli d4,d2,e
  ry f,0.007569891767651322
  cx e,f
  ry f,-0.007569891767651322
  cx d1,e
  cx d1,f
  ry f,0.007569891767651322
  cx e,f
  cx d1,e
  cx d1,f
  ry f,-0.007569891767651322
  cx d1,f
  toffoli d4,d2,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.0234375, controlled by [0, 3, 5]
  toffoli d4,d1,e
  ry f,0.0009462364709564153
  cx e,f
  ry f,-0.0009462364709564153
  cx c5,e
  cx c5,f
  ry f,0.0009462364709564153
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.0009462364709564153
  cx c5,f
  toffoli d4,d1,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.375, controlled by [0, 1, 3]
  toffoli d4,d3,e
  ry f,0.015139783535302645
  cx e,f
  ry f,-0.015139783535302645
  cx d1,e
  cx d1,f
  ry f,0.015139783535302645
  cx e,f
  cx d1,e
  cx d1,f
  ry f,-0.015139783535302645
  cx d1,f
  toffoli d4,d3,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.09375, controlled by [1, 2, 3]
  toffoli d3,d2,e
  ry f,0.003784945883825661
  cx e,f
  ry f,-0.003784945883825661
  cx d1,e
  cx d1,f
  ry f,0.003784945883825661
  cx e,f
  cx d1,e
  cx d1,f
  ry f,-0.003784945883825661
  cx d1,f
  toffoli d3,d2,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.0234375, controlled by [1, 3, 4]
  toffoli d3,d1,e
  ry f,0.0009462364709564153
  cx e,f
  ry f,-0.0009462364709564153
  cx d0,e
  cx d0,f
  ry f,0.0009462364709564153
  cx e,f
  cx d0,e
  cx d0,f
  ry f,-0.0009462364709564153
  cx d0,f
  toffoli d3,d1,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.01171875, controlled by [1, 3, 5]
  toffoli d3,d1,e
  ry f,0.00047311823547820764
  cx e,f
  ry f,-0.00047311823547820764
  cx c5,e
  cx c5,f
  ry f,0.00047311823547820764
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.00047311823547820764
  cx c5,f
  toffoli d3,d1,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.005859375, controlled by [1, 4, 5]
  toffoli d3,d0,e
  ry f,0.00023655911773910382
  cx e,f
  ry f,-0.00023655911773910382
  cx c5,e
  cx c5,f
  ry f,0.00023655911773910382
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.00023655911773910382
  cx c5,f
  toffoli d3,d0,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.1875, controlled by [0, 1, 4]
  toffoli d4,d3,e
  ry f,0.007569891767651322
  cx e,f
  ry f,-0.007569891767651322
  cx d0,e
  cx d0,f
  ry f,0.007569891767651322
  cx e,f
  cx d0,e
  cx d0,f
  ry f,-0.007569891767651322
  cx d0,f
  toffoli d4,d3,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.046875, controlled by [0, 2, 5]
  toffoli d4,d2,e
  ry f,0.0018924729419128306
  cx e,f
  ry f,-0.0018924729419128306
  cx c5,e
  cx c5,f
  ry f,0.0018924729419128306
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.0018924729419128306
  cx c5,f
  toffoli d4,d2,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.046875, controlled by [1, 2, 4]
  toffoli d3,d2,e
  ry f,0.0018924729419128306
  cx e,f
  ry f,-0.0018924729419128306
  cx d0,e
  cx d0,f
  ry f,0.0018924729419128306
  cx e,f
  cx d0,e
  cx d0,f
  ry f,-0.0018924729419128306
  cx d0,f
  toffoli d3,d2,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.0029296875, controlled by [2, 4, 5]
  toffoli d2,d0,e
  ry f,0.00011827955886955191
  cx e,f
  ry f,-0.00011827955886955191
  cx c5,e
  cx c5,f
  ry f,0.00011827955886955191
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.00011827955886955191
  cx c5,f
  toffoli d2,d0,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.01171875, controlled by [2, 3, 4]
  toffoli d2,d1,e
  ry f,0.00047311823547820764
  cx e,f
  ry f,-0.00047311823547820764
  cx d0,e
  cx d0,f
  ry f,0.00047311823547820764
  cx e,f
  cx d0,e
  cx d0,f
  ry f,-0.00047311823547820764
  cx d0,f
  toffoli d2,d1,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.09375, controlled by [0, 2, 4]
  toffoli d4,d2,e
  ry f,0.003784945883825661
  cx e,f
  ry f,-0.003784945883825661
  cx d0,e
  cx d0,f
  ry f,0.003784945883825661
  cx e,f
  cx d0,e
  cx d0,f
  ry f,-0.003784945883825661
  cx d0,f
  toffoli d4,d2,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.0234375, controlled by [1, 2, 5]
  toffoli d3,d2,e
  ry f,0.0009462364709564153
  cx e,f
  ry f,-0.0009462364709564153
  cx c5,e
  cx c5,f
  ry f,0.0009462364709564153
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.0009462364709564153
  cx c5,f
  toffoli d3,d2,e
  
  # Rotation c*theta with c = 0.08074551218828077 and theta = 0.09375, controlled by [0, 1, 5]
  toffoli d4,d3,e
  ry f,0.003784945883825661
  cx e,f
  ry f,-0.003784945883825661
  cx c5,e
  cx c5,f
  ry f,0.003784945883825661
  cx e,f
  cx c5,e
  cx c5,f
  ry f,-0.003784945883825661
  cx c5,f
  toffoli d4,d3,e
  

.thapliyal_division
  ##### un-Thapliyal division with Q = [c0, ..., c5], R = [d0, ..., d4] and D = [b0, ..., b4]
  
  ### Un-iteration 6, where Y = [c0, ..., c4] and R = c5
  
  # x on Z
  x c5
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x c0
  x c1
  x c2
  x c3
  x c4
  cx b1,c1
  cx b2,c2
  cx b3,c3
  cx b4,c4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c0,b1
  toffoli b1,c1,b2
  toffoli b2,c2,b3
  toffoli b3,c3,b4
  toffoli c5,b4,c4
  toffoli b3,c3,b4
  toffoli c5,b3,c3
  toffoli b2,c2,b3
  toffoli c5,b2,c2
  toffoli b1,c1,b2
  toffoli c5,b1,c1
  toffoli b0,c0,b1
  toffoli c5,b0,c0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c4
  cx b3,c3
  cx b2,c2
  cx b1,c1
  x c0
  x c1
  x c2
  x c3
  x c4
  
  # un-SUB: Y - D, on Y
  cx b1,c1
  cx b2,c2
  cx b3,c3
  cx b4,c4
  cx b4,c5
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c0,b1
  toffoli b1,c1,b2
  toffoli b2,c2,b3
  toffoli b3,c3,b4
  toffoli b4,c4,c5
  cx b4,c4
  toffoli b3,c3,b4
  cx b3,c3
  toffoli b2,c2,b3
  cx b2,c2
  toffoli b1,c1,b2
  cx b1,c1
  toffoli b0,c0,b1
  cx b0,c0
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c4
  cx b3,c3
  cx b2,c2
  cx b1,c1
  
  ### Un-iteration 5, where Y = [c1, ..., c5] and R = d0
  
  # x on Z
  x d0
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x c1
  x c2
  x c3
  x c4
  x c5
  cx b1,c2
  cx b2,c3
  cx b3,c4
  cx b4,c5
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c1,b1
  toffoli b1,c2,b2
  toffoli b2,c3,b3
  toffoli b3,c4,b4
  toffoli d0,b4,c5
  toffoli b3,c4,b4
  toffoli d0,b3,c4
  toffoli b2,c3,b3
  toffoli d0,b2,c3
  toffoli b1,c2,b2
  toffoli d0,b1,c2
  toffoli b0,c1,b1
  toffoli d0,b0,c1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c5
  cx b3,c4
  cx b2,c3
  cx b1,c2
  x c1
  x c2
  x c3
  x c4
  x c5
  
  # un-SUB: Y - D, on Y
  cx b1,c2
  cx b2,c3
  cx b3,c4
  cx b4,c5
  cx b4,d0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c1,b1
  toffoli b1,c2,b2
  toffoli b2,c3,b3
  toffoli b3,c4,b4
  toffoli b4,c5,d0
  cx b4,c5
  toffoli b3,c4,b4
  cx b3,c4
  toffoli b2,c3,b3
  cx b2,c3
  toffoli b1,c2,b2
  cx b1,c2
  toffoli b0,c1,b1
  cx b0,c1
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,c5
  cx b3,c4
  cx b2,c3
  cx b1,c2
  
  ### Un-iteration 4, where Y = [c2, ..., d0] and R = d1
  
  # x on Z
  x d1
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x c2
  x c3
  x c4
  x c5
  x d0
  cx b1,c3
  cx b2,c4
  cx b3,c5
  cx b4,d0
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c2,b1
  toffoli b1,c3,b2
  toffoli b2,c4,b3
  toffoli b3,c5,b4
  toffoli d1,b4,d0
  toffoli b3,c5,b4
  toffoli d1,b3,c5
  toffoli b2,c4,b3
  toffoli d1,b2,c4
  toffoli b1,c3,b2
  toffoli d1,b1,c3
  toffoli b0,c2,b1
  toffoli d1,b0,c2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d0
  cx b3,c5
  cx b2,c4
  cx b1,c3
  x c2
  x c3
  x c4
  x c5
  x d0
  
  # un-SUB: Y - D, on Y
  cx b1,c3
  cx b2,c4
  cx b3,c5
  cx b4,d0
  cx b4,d1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c2,b1
  toffoli b1,c3,b2
  toffoli b2,c4,b3
  toffoli b3,c5,b4
  toffoli b4,d0,d1
  cx b4,d0
  toffoli b3,c5,b4
  cx b3,c5
  toffoli b2,c4,b3
  cx b2,c4
  toffoli b1,c3,b2
  cx b1,c3
  toffoli b0,c2,b1
  cx b0,c2
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d0
  cx b3,c5
  cx b2,c4
  cx b1,c3
  
  ### Un-iteration 3, where Y = [c3, ..., d1] and R = d2
  
  # x on Z
  x d2
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x c3
  x c4
  x c5
  x d0
  x d1
  cx b1,c4
  cx b2,c5
  cx b3,d0
  cx b4,d1
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c3,b1
  toffoli b1,c4,b2
  toffoli b2,c5,b3
  toffoli b3,d0,b4
  toffoli d2,b4,d1
  toffoli b3,d0,b4
  toffoli d2,b3,d0
  toffoli b2,c5,b3
  toffoli d2,b2,c5
  toffoli b1,c4,b2
  toffoli d2,b1,c4
  toffoli b0,c3,b1
  toffoli d2,b0,c3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d1
  cx b3,d0
  cx b2,c5
  cx b1,c4
  x c3
  x c4
  x c5
  x d0
  x d1
  
  # un-SUB: Y - D, on Y
  cx b1,c4
  cx b2,c5
  cx b3,d0
  cx b4,d1
  cx b4,d2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c3,b1
  toffoli b1,c4,b2
  toffoli b2,c5,b3
  toffoli b3,d0,b4
  toffoli b4,d1,d2
  cx b4,d1
  toffoli b3,d0,b4
  cx b3,d0
  toffoli b2,c5,b3
  cx b2,c5
  toffoli b1,c4,b2
  cx b1,c4
  toffoli b0,c3,b1
  cx b0,c3
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d1
  cx b3,d0
  cx b2,c5
  cx b1,c4
  
  ### Un-iteration 2, where Y = [c4, ..., d2] and R = d3
  
  # x on Z
  x d3
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x c4
  x c5
  x d0
  x d1
  x d2
  cx b1,c5
  cx b2,d0
  cx b3,d1
  cx b4,d2
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c4,b1
  toffoli b1,c5,b2
  toffoli b2,d0,b3
  toffoli b3,d1,b4
  toffoli d3,b4,d2
  toffoli b3,d1,b4
  toffoli d3,b3,d1
  toffoli b2,d0,b3
  toffoli d3,b2,d0
  toffoli b1,c5,b2
  toffoli d3,b1,c5
  toffoli b0,c4,b1
  toffoli d3,b0,c4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d2
  cx b3,d1
  cx b2,d0
  cx b1,c5
  x c4
  x c5
  x d0
  x d1
  x d2
  
  # un-SUB: Y - D, on Y
  cx b1,c5
  cx b2,d0
  cx b3,d1
  cx b4,d2
  cx b4,d3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c4,b1
  toffoli b1,c5,b2
  toffoli b2,d0,b3
  toffoli b3,d1,b4
  toffoli b4,d2,d3
  cx b4,d2
  toffoli b3,d1,b4
  cx b3,d1
  toffoli b2,d0,b3
  cx b2,d0
  toffoli b1,c5,b2
  cx b1,c5
  toffoli b0,c4,b1
  cx b0,c4
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d2
  cx b3,d1
  cx b2,d0
  cx b1,c5
  
  ### Un-iteration 1, where Y = [c5, ..., d3] and R = d4
  
  # x on Z
  x d4
  
  # un-cADD: Y (+ D), on Y, ctrl'd by Z
  x c5
  x d0
  x d1
  x d2
  x d3
  cx b1,d0
  cx b2,d1
  cx b3,d2
  cx b4,d3
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c5,b1
  toffoli b1,d0,b2
  toffoli b2,d1,b3
  toffoli b3,d2,b4
  toffoli d4,b4,d3
  toffoli b3,d2,b4
  toffoli d4,b3,d2
  toffoli b2,d1,b3
  toffoli d4,b2,d1
  toffoli b1,d0,b2
  toffoli d4,b1,d0
  toffoli b0,c5,b1
  toffoli d4,b0,c5
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d3
  cx b3,d2
  cx b2,d1
  cx b1,d0
  x c5
  x d0
  x d1
  x d2
  x d3
  
  # un-SUB: Y - D, on Y
  cx b1,d0
  cx b2,d1
  cx b3,d2
  cx b4,d3
  cx b4,d4
  cx b3,b4
  cx b2,b3
  cx b1,b2
  toffoli b0,c5,b1
  toffoli b1,d0,b2
  toffoli b2,d1,b3
  toffoli b3,d2,b4
  toffoli b4,d3,d4
  cx b4,d3
  toffoli b3,d2,b4
  cx b3,d2
  toffoli b2,d1,b3
  cx b2,d1
  toffoli b1,d0,b2
  cx b1,d0
  toffoli b0,c5,b1
  cx b0,c5
  cx b1,b2
  cx b2,b3
  cx b3,b4
  cx b4,d3
  cx b3,d2
  cx b2,d1
  cx b1,d0

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
  x c5

.result
  display
  measure
  display