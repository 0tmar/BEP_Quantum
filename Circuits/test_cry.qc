#function: cRy test

qubits 2

.test

  # init
  x q0
  x q1
  
  # cRy(a,b,-pi/2)
  ry q1,-0.19634954084936207
  cx q0,q1
  ry q1,0.39269908169872414
  cx q0,q1
  ry q1,-0.19634954084936207
  display
  
  # Ry(b,pi/2)
  ry q1,0.7853981633974483
  display