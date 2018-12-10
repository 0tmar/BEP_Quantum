#function: ccZ test

qubits 3

.test

  # init
  h q0
  h q1
  h q2
  display
  
  # ccZ(a,b,c)
  tdag q1
  cx q1,q2
  t q2
  cx q1,q2
  cx q0,q1
  t q1
  cx q1,q2
  tdag q2
  cx q1,q2
  cx q0,q1
  tdag q0
  tdag q2
  cx q0,q2
  t q2
  cx q0,q2
  display
  
  # ccZ(a,b,c)
  tdag q1
  cx q1,q2
  t q2
  cx q1,q2
  cx q0,q1
  t q1
  cx q1,q2
  tdag q2
  cx q1,q2
  cx q0,q1
  tdag q0
  tdag q2
  cx q0,q2
  t q2
  cx q0,q2
  display

  # init
  h q0
  h q1
  h q2
  display
  measure
  display