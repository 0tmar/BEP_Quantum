#function: Cao2012 Experiment

qubits 7

.init
  h q1
  h q2
  h q3
  h q4
  h q5
  h q6
  display

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/16) on q5 and q6, controlled by q1
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q1,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q1,q5
  tdag q1
  tdag q6
  cx q1,q6
  t q6
  cx q1,q6
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q1,q6
  rx q6,0.098
  cz q1,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q1
  x q1
  tdag q1
  x q1
  h q6
  tdag q6
  tdag q6
  cx q1,q6
  h q6
  t q6
  cx q1,q6
  t q6
  h q6
  s q6
  x q6
  tdag q1
  tdag q6
  cx q1,q6
  t q6
  cx q1,q6
  t q1
  x q1
  tdag q1
  x q1
  h q6
  tdag q6
  tdag q6
  cx q1,q6
  h q6
  t q6
  cx q1,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q1,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q1,q5
  rx q5,0.491
  cz q1,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q1,1.883
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q1,q5
  rx q5,0.2945
  cz q1,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q1,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q1,q5
  tdag q1
  tdag q6
  cx q1,q6
  t q6
  cx q1,q6

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/8) on q5 and q6, controlled by q2
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q2,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q2,q5
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6
  
  # cRx(a,c,-0.392)
  rx q6,-0.098
  cz q2,q6
  rx q6,0.196
  cz q2,q6
  rx q6,-0.098
  
  # cX(a,c)
  cx q2,q6
  
  # Rz(a,0.75)
  rz q2,0.75
  
  # cRx(a,b,-1.964)
  rx q5,-0.491
  cz q2,q5
  rx q5,0.982
  cz q2,q5
  rx q5,-0.491
  
  # Rz(a,3.766)
  rz q2,3.766
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,-1.178)
  rx q5,-0.2945
  cz q2,q5
  rx q5,0.589
  cz q2,q5
  rx q5,-0.2945
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q2,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q2,q5
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/4) on q5 and q6, controlled by q3
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q3,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q3,q5
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  
  # cRx(a,c,-0.784)
  rx q6,-0.196
  cz q3,q6
  rx q6,0.392
  cz q3,q6
  rx q6,-0.196
  
  # Rz(a,1.5)
  rz q3,1.5
  
  # cRx(a,b,-3.928)
  rx q5,-0.982
  cz q3,q5
  rx q5,1.964
  cz q3,q5
  rx q5,-0.982
  
  # Rz(a,7.532)
  rz q3,7.532
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,-2.356)
  rx q5,-0.589
  cz q3,q5
  rx q5,1.178
  cz q3,q5
  rx q5,-0.589
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q3,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q3,q5
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/2) on q5 and q6, controlled by q4
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q4,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q4,q5
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  
  # cRx(a,c,-1.568)
  rx q6,-0.392
  cz q4,q6
  rx q6,0.784
  cz q4,q6
  rx q6,-0.392
  
  # Rz(a,3.0)
  rz q4,3.0
  
  # cRx(a,b,-7.856)
  rx q5,-1.964
  cz q4,q5
  rx q5,3.928
  cz q4,q5
  rx q5,-1.964
  
  # Rz(a,15.064)
  rz q4,15.064
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-4.712)
  rx q5,-1.178
  cz q4,q5
  rx q5,2.356
  cz q4,q5
  rx q5,-1.178
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q4,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q4,q5
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6

.iqft
  display
  
  h q4
  cr q4,q3
  cz q4,q3
  h q3
  cr q4,q2
  swap q4,q3
  cr q3,q2
  swap q4,q3
  cz q4,q2
  cr q3,q2
  cz q3,q2
  h q2
  cr q4,q1
  swap q4,q3
  cr q3,q1
  swap q4,q3
  swap q4,q2
  cr q2,q1
  swap q4,q2
  cz q4,q1
  cr q3,q1
  swap q3,q2
  cr q2,q1
  swap q3,q2
  cz q3,q1
  cr q2,q1
  cz q2,q1
  h q1
  
  display

.reverse
  swap q1,q4
  swap q2,q3

.cRy
  # performs Ry(2^(3)*pi/2^(6-1)) = Ry(pi/4) on q0, controlled by q1
  ry q0,0.19634954084936207
  cx q1,q0
  ry q0,-0.39269908169872414
  cx q1,q0
  ry q0,0.19634954084936207

.cRy
  # performs Ry(2^(2)*pi/2^(6-1)) = Ry(pi/8) on q0, controlled by q2
  ry q0,0.09817477042468103
  cx q2,q0
  ry q0,-0.19634954084936207
  cx q2,q0
  ry q0,0.09817477042468103

.cRy
  # performs Ry(2^(1)*pi/2^(6-1)) = Ry(pi/16) on q0, controlled by q3
  ry q0,0.04908738521234052
  cx q3,q0
  ry q0,-0.09817477042468103
  cx q3,q0
  ry q0,0.04908738521234052

.cRy
  # performs Ry(2^(0)*pi/2^(6-1)) = Ry(pi/32) on q0, controlled by q4
  ry q0,0.02454369260617026
  cx q4,q0
  ry q0,-0.04908738521234052
  cx q4,q0
  ry q0,0.02454369260617026

.reverse
  swap q1,q4
  swap q2,q3

.qft
  h q1
  cr q2,q1
  cr q3,q1
  cr q4,q1
  h q2
  cr q3,q2
  cr q4,q2
  h q3
  cr q4,q3
  h q4

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/2) on q5 and q6, controlled by q4
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q4,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q4,q5
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,4.712)
  rx q5,1.178
  cz q4,q5
  rx q5,-2.356
  cz q4,q5
  rx q5,1.178
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-15.064)
  rz q4,-15.064
  
  # cRx(a,b,7.856)
  rx q5,1.964
  cz q4,q5
  rx q5,-3.928
  cz q4,q5
  rx q5,1.964
  
  # Rz(a,-3.0)
  rz q4,-3.0
  
  # cRx(a,c,1.568)
  rx q6,0.392
  cz q4,q6
  rx q6,-0.784
  cz q4,q6
  rx q6,0.392
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q4,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q4,q5
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/4) on q5 and q6, controlled by q3
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q3,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q3,q5
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,2.356)
  rx q5,0.589
  cz q3,q5
  rx q5,-1.178
  cz q3,q5
  rx q5,0.589
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # Rz(a,-7.532)
  rz q3,-7.532
  
  # cRx(a,b,3.928)
  rx q5,0.982
  cz q3,q5
  rx q5,-1.964
  cz q3,q5
  rx q5,0.982
  
  # Rz(a,-1.5)
  rz q3,-1.5
  
  # cRx(a,c,0.784)
  rx q6,0.196
  cz q3,q6
  rx q6,-0.392
  cz q3,q6
  rx q6,0.196
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q3,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q3,q5
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/8) on q5 and q6, controlled by q2
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q2,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q2,q5
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,1.178)
  rx q5,0.2945
  cz q2,q5
  rx q5,-0.589
  cz q2,q5
  rx q5,0.2945
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # Rz(a,-3.766)
  rz q2,-3.766
  
  # cRx(a,b,1.964)
  rx q5,0.491
  cz q2,q5
  rx q5,-0.982
  cz q2,q5
  rx q5,0.491
  
  # Rz(a,-0.75)
  rz q2,-0.75
  
  # cX(a,c)
  cx q2,q6
  
  # cRx(a,c,0.392)
  rx q6,0.098
  cz q2,q6
  rx q6,-0.196
  cz q2,q6
  rx q6,0.098
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q2,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q2,q5
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/16) on q5 and q6, controlled by q1
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q1,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q1,q5
  tdag q1
  tdag q6
  cx q1,q6
  t q6
  cx q1,q6
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q1,q5
  rx q5,-0.2945
  cz q1,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # Rz(a,-1.883)
  rz q1,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q1,q5
  rx q5,-0.491
  cz q1,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q1,-0.375
  
  # cV(a,c)
  t q1
  x q1
  tdag q1
  x q1
  h q6
  tdag q6
  tdag q6
  cx q1,q6
  h q6
  t q6
  cx q1,q6
  t q6
  h q6
  s q6
  x q6
  t q1
  t q6
  cx q1,q6
  tdag q6
  cx q1,q6
  t q1
  x q1
  tdag q1
  x q1
  h q6
  tdag q6
  tdag q6
  cx q1,q6
  h q6
  t q6
  cx q1,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q1,q6
  rx q6,-0.098
  cz q1,q6
  rx q6,0.049
  
  # ccZ(a,b,c)
  tdag q5
  cx q5,q6
  t q6
  cx q5,q6
  cx q1,q5
  t q5
  cx q5,q6
  tdag q6
  cx q5,q6
  cx q1,q5
  tdag q1
  tdag q6
  cx q1,q6
  t q6
  cx q1,q6

.uninit
  h q1
  h q2
  h q3
  h q4

.result
  display
  measure
  display