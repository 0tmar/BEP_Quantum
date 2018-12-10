#function: cH test

qubits 2

.test

  # init
  x q0
  x q1
  display
  
  # cH(a,b)
  t q0
  x q0
  tdag q0
  x q0
  h q1
  tdag q1
  tdag q1
  cx q0,q1
  h q1
  t q1
  cx q0,q1
  t q1
  h q1
  s q1
  x q1
  display
  
  # cH(a,b)
  t q0
  x q0
  tdag q0
  x q0
  h q1
  tdag q1
  tdag q1
  cx q0,q1
  h q1
  t q1
  cx q0,q1
  t q1
  h q1
  s q1
  x q1
  display

  # meas
  measure
  display