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
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/16) on q5 and q6, controlled by q1
  
  # cX(a,b)
  cx q1,q5
  
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
  
  # cX(a,b)
  cx q1,q5
  
  # cZ(a,c)
  cz q1,q6
  
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

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/8) on q5 and q6, controlled by q2
  
  # cX(a,b)
  cx q2,q5
  
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
  
  # cX(a,b)
  cx q2,q5
  
  # cZ(a,c)
  cz q2,q6
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q2,q5
  rx q5,-0.2945
  cz q2,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # Rz(a,-1.883)
  rz q2,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q2,q5
  rx q5,-0.491
  cz q2,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q2,-0.375
  
  # cV(a,c)
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  t q2
  t q6
  cx q2,q6
  tdag q6
  cx q2,q6
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q2,q6
  rx q6,-0.098
  cz q2,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q2,q5
  
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
  
  # cX(a,b)
  cx q2,q5
  
  # cZ(a,c)
  cz q2,q6
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q2,q5
  rx q5,-0.2945
  cz q2,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # Rz(a,-1.883)
  rz q2,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q2,q5
  rx q5,-0.491
  cz q2,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q2,-0.375
  
  # cV(a,c)
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  t q2
  t q6
  cx q2,q6
  tdag q6
  cx q2,q6
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q2,q6
  rx q6,-0.098
  cz q2,q6
  rx q6,0.049
  
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
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/4) on q5 and q6, controlled by q3
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5
  
  # cZ(a,c)
  cz q3,q6
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q3,q5
  rx q5,-0.2945
  cz q3,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # Rz(a,-1.883)
  rz q3,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q3,q5
  rx q5,-0.491
  cz q3,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q3,-0.375
  
  # cV(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q3,q6
  rx q6,-0.098
  cz q3,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5
  
  # cZ(a,c)
  cz q3,q6
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q3,q5
  rx q5,-0.2945
  cz q3,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # Rz(a,-1.883)
  rz q3,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q3,q5
  rx q5,-0.491
  cz q3,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q3,-0.375
  
  # cV(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q3,q6
  rx q6,-0.098
  cz q3,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5
  
  # cZ(a,c)
  cz q3,q6
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q3,q5
  rx q5,-0.2945
  cz q3,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # Rz(a,-1.883)
  rz q3,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q3,q5
  rx q5,-0.491
  cz q3,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q3,-0.375
  
  # cV(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q3,q6
  rx q6,-0.098
  cz q3,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5
  
  # cZ(a,c)
  cz q3,q6
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q3,q5
  rx q5,-0.2945
  cz q3,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # Rz(a,-1.883)
  rz q3,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q3,q5
  rx q5,-0.491
  cz q3,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q3,-0.375
  
  # cV(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q3,q6
  rx q6,-0.098
  cz q3,q6
  rx q6,0.049
  
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
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/2) on q5 and q6, controlled by q4
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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
  
  #####
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  # cZ(a,c)
  cz q4,q6
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.589)
  rx q5,0.14725
  cz q4,q5
  rx q5,-0.2945
  cz q4,q5
  rx q5,0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # Rz(a,-1.883)
  rz q4,-1.883
  
  # cRx(a,b,0.982)
  rx q5,0.2455
  cz q4,q5
  rx q5,-0.491
  cz q4,q5
  rx q5,0.2455
  
  # Rz(a,-0.375)
  rz q4,-0.375
  
  # cV(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # cRx(a,c,0.196)
  rx q6,0.049
  cz q4,q6
  rx q6,-0.098
  cz q4,q6
  rx q6,0.049
  
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

.swap
  swap q2,q4

.cRy
  # performs Ry(-2^(3)*pi/2^(5-1)) = Ry(-pi/2) on q0, controlled by q1
  ry q0,-0.39269908169872414
  cx q1,q0
  ry q0,0.7853981633974483
  cx q1,q0
  ry q0,-0.39269908169872414

.cRy
  # performs Ry(-2^(2)*pi/2^(5-1)) = Ry(-pi/4) on q0, controlled by q2
  ry q0,-0.19634954084936207
  cx q2,q0
  ry q0,0.39269908169872414
  cx q2,q0
  ry q0,-0.19634954084936207

.cRy
  # performs Ry(-2^(1)*pi/2^(5-1)) = Ry(-pi/8) on q0, controlled by q3
  ry q0,-0.09817477042468103
  cx q3,q0
  ry q0,0.19634954084936207
  cx q3,q0
  ry q0,-0.09817477042468103

.cRy
  # performs Ry(-2^(0)*pi/2^(5-1)) = Ry(-pi/16) on q0, controlled by q4
  ry q0,-0.04908738521234052
  cx q4,q0
  ry q0,0.09817477042468103
  cx q4,q0
  ry q0,-0.04908738521234052

.swap
  swap q2,q4

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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q4,q6
  rx q6,0.098
  cz q4,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  cz q4,q6
  t q4
  t q6
  cx q4,q6
  tdag q6
  cx q4,q6
  t q4
  x q4
  tdag q4
  x q4
  h q6
  tdag q6
  tdag q6
  cx q4,q6
  h q6
  t q6
  cx q4,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q4,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q4,q5
  rx q5,0.491
  cz q4,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q4,1.883
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q4,q5
  rx q5,0.2945
  cz q4,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
  # cX(a,b)
  cx q4,q5
  
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
  
  # cX(a,b)
  cx q4,q5

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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q3,q6
  rx q6,0.098
  cz q3,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  cz q3,q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q3,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q3,q5
  rx q5,0.491
  cz q3,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q3,1.883
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q3,q5
  rx q5,0.2945
  cz q3,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q3,q6
  rx q6,0.098
  cz q3,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  cz q3,q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q3,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q3,q5
  rx q5,0.491
  cz q3,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q3,1.883
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q3,q5
  rx q5,0.2945
  cz q3,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q3,q6
  rx q6,0.098
  cz q3,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  cz q3,q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q3,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q3,q5
  rx q5,0.491
  cz q3,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q3,1.883
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q3,q5
  rx q5,0.2945
  cz q3,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q3,q6
  rx q6,0.098
  cz q3,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  cz q3,q6
  t q3
  t q6
  cx q3,q6
  tdag q6
  cx q3,q6
  t q3
  x q3
  tdag q3
  x q3
  h q6
  tdag q6
  tdag q6
  cx q3,q6
  h q6
  t q6
  cx q3,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q3,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q3,q5
  rx q5,0.491
  cz q3,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q3,1.883
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q3,q5
  rx q5,0.2945
  cz q3,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
  # cX(a,b)
  cx q3,q5
  
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
  
  # cX(a,b)
  cx q3,q5

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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q2,q6
  rx q6,0.098
  cz q2,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  cz q2,q6
  t q2
  t q6
  cx q2,q6
  tdag q6
  cx q2,q6
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q2,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q2,q5
  rx q5,0.491
  cz q2,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q2,1.883
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q2,q5
  rx q5,0.2945
  cz q2,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cZ(a,c)
  cz q2,q6
  
  # cX(a,b)
  cx q2,q5
  
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
  
  # cX(a,b)
  cx q2,q5
  
  #####
  
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
  
  # cRx(a,c,-0.196)
  rx q6,-0.049
  cz q2,q6
  rx q6,0.098
  cz q2,q6
  rx q6,-0.049
  
  # cVdag(a,c)
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  cz q2,q6
  t q2
  t q6
  cx q2,q6
  tdag q6
  cx q2,q6
  t q2
  x q2
  tdag q2
  x q2
  h q6
  tdag q6
  tdag q6
  cx q2,q6
  h q6
  t q6
  cx q2,q6
  t q6
  h q6
  s q6
  x q6
  
  # Rz(a,0.375)
  rz q2,0.375
  
  # cRx(a,b,-0.982)
  rx q5,-0.2455
  cz q2,q5
  rx q5,0.491
  cz q2,q5
  rx q5,-0.2455
  
  # Rz(a,1.883)
  rz q2,1.883
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,-0.589)
  rx q5,-0.14725
  cz q2,q5
  rx q5,0.2945
  cz q2,q5
  rx q5,-0.14725
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cZ(a,c)
  cz q2,q6
  
  # cX(a,b)
  cx q2,q5
  
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
  
  # cX(a,b)
  cx q2,q5

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
  cz q1,q6
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
  
  # cZ(a,c)
  cz q1,q6
  
  # cX(a,b)
  cx q1,q5
  
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
  
  # cX(a,b)
  cx q1,q5

.uninit
  h q1
  h q2
  h q3
  h q4

.result
  display
  measure
  display