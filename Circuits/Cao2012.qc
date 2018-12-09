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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q1,q6
  rx q6,-0.1
  cz q1,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q1
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
  tdag q1
  tdag q6
  cx q1,q6
  t q6
  cx q1,q6
  s q1
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
  
  # Rz(a,0.38)
  rz q1,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q1,q6
  rx q6,-0.49
  cz q1,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q1,1.88
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q1,q6
  rx q6,-0.295
  cz q1,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # cZ(a,c)
  cz q1,q6

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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q2,q6
  rx q6,-0.1
  cz q2,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q2
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
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6
  s q2
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
  
  # Rz(a,0.38)
  rz q2,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q2,q6
  rx q6,-0.49
  cz q2,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q2,1.88
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q2,q6
  rx q6,-0.295
  cz q2,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cZ(a,c)
  cz q2,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q2,q6
  rx q6,-0.1
  cz q2,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q2
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
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6
  s q2
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
  
  # Rz(a,0.38)
  rz q2,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q2,q6
  rx q6,-0.49
  cz q2,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q2,1.88
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q2,q6
  rx q6,-0.295
  cz q2,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cZ(a,c)
  cz q2,q6

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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6

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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6

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
  # performs Ry(2^(3)*pi/2^(5-1)) = Ry(pi/2^2) on q0, controlled by q1
  ry q0,0.39269908169872414
  cx q1,q0
  ry q0,-0.7853981633974483
  cx q1,q0
  ry q0,0.39269908169872414

.cRy
  # performs Ry(2^(2)*pi/2^(5-1)) = Ry(pi/2^4) on q0, controlled by q2
  ry q0,0.19634954084936207
  cx q2,q0
  ry q0,-0.39269908169872414
  cx q2,q0
  ry q0,0.19634954084936207

.cRy
  # performs Ry(2^(1)*pi/2^(5-1)) = Ry(pi/2^8) on q0, controlled by q3
  ry q0,0.09817477042468103
  cx q3,q0
  ry q0,-0.19634954084936207
  cx q3,q0
  ry q0,0.09817477042468103

.cRy
  # performs Ry(2^(0)*pi/2^(5-1)) = Ry(pi/2^16) on q0, controlled by q4
  ry q0,0.04908738521234052
  cx q4,q0
  ry q0,-0.09817477042468103
  cx q4,q0
  ry q0,0.04908738521234052

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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q4,q6
  rx q6,-0.1
  cz q4,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q4
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
  tdag q4
  tdag q6
  cx q4,q6
  t q6
  cx q4,q6
  s q4
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
  
  # Rz(a,0.38)
  rz q4,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q4,q6
  rx q6,-0.49
  cz q4,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q4,1.88
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q4,q6
  rx q6,-0.295
  cz q4,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q4,q5,q6
  
  # cZ(a,c)
  cz q4,q6

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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q3,q6
  rx q6,-0.1
  cz q3,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q3
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
  tdag q3
  tdag q6
  cx q3,q6
  t q6
  cx q3,q6
  s q3
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
  
  # Rz(a,0.38)
  rz q3,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q3,q6
  rx q6,-0.49
  cz q3,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q3,1.88
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q3,q6
  rx q6,-0.295
  cz q3,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q3,q5,q6
  
  # cZ(a,c)
  cz q3,q6

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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q2,q6
  rx q6,-0.1
  cz q2,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q2
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
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6
  s q2
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
  
  # Rz(a,0.38)
  rz q2,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q2,q6
  rx q6,-0.49
  cz q2,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q2,1.88
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q2,q6
  rx q6,-0.295
  cz q2,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cZ(a,c)
  cz q2,q6
  
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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q2,q6
  rx q6,-0.1
  cz q2,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q2
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
  tdag q2
  tdag q6
  cx q2,q6
  t q6
  cx q2,q6
  s q2
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
  
  # Rz(a,0.38)
  rz q2,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q2,q6
  rx q6,-0.49
  cz q2,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q2,1.88
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q2,q6
  rx q6,-0.295
  cz q2,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q2,q5,q6
  
  # cZ(a,c)
  cz q2,q6

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
  
  # cRx(a,c,0.2)
  rx q6,0.05
  cz q1,q6
  rx q6,-0.1
  cz q1,q6
  rx q6,0.05
  
  # cVdag(a,c)
  s q1
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
  tdag q1
  tdag q6
  cx q1,q6
  t q6
  cx q1,q6
  s q1
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
  
  # Rz(a,0.38)
  rz q1,0.38
  
  # cRx(a,b,0.98)
  rx q6,0.245
  cz q1,q6
  rx q6,-0.49
  cz q1,q6
  rx q6,0.245
  
  # Rz(a,1.88)
  rz q1,1.88
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # cRx(a,b,0.59)
  rx q6,0.1475
  cz q1,q6
  rx q6,-0.295
  cz q1,q6
  rx q6,0.1475
  
  # Toffoli(a,b,c)
  toffoli q1,q5,q6
  
  # cZ(a,c)
  cz q1,q6

.uninit
  h q1
  h q2
  h q3
  h q4

.result
  display
  measure
  display