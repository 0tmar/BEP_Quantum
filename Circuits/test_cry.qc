#function: cRy test

qubits 2

.test

  # init
  x q0
  x q1
  display
  
  # cRy(a,b,0.2)
  ry q1,0.05
  cx q0,q1
  ry q1,-0.1
  cx q0,q1
  ry q1,0.05
  display
  
  ## cRy(a,b,-0.2)
  #ry q1,-0.05
  #cx q0,q1
  #ry q1,0.1
  #cx q0,q1
  #ry q1,-0.05
  #display
  
  # Ry(b,-0.2)
  ry q1,-0.2
  display

  # init
  measure
  display