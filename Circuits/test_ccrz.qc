#function: ccRz test

qubits 3

.test
  
  # init
  x q0
  x q1
  # x q2
  display
  
  # ccRz(a,b,c,4.0)
  rz q2,1.0
  cx q1,q2
  rz q2,-1.0
  cx q0,q1
  cx q0,q2
  rz q2,1.0
  cx q1,q2
  cx q0,q1
  cx q0,q2
  rz q2,-1.0
  cx q0,q2
  display
  
  # meas
  measure
  display