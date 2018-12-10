#function: cS-cSdag test

qubits 2

.test

  # init
  x q0
  x q1
  display
  
  # cS(a,b)
  t q0
  t q1
  cx q0,q1
  tdag q1
  cx q0,q1
  display
  
  # cSdag(a,b)
  tdag q0
  tdag q1
  cx q0,q1
  t q1
  cx q0,q1
  display

  # meas
  measure
  display