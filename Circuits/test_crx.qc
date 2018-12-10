#function: cRx test

qubits 2

.test

  # init
  x q0
  x q1
  display
  
  # cRx(a,b,0.2)
  rx q1,0.05
  cz q0,q1
  rx q1,-0.1
  cz q0,q1
  rx q1,0.05
  display
  
  ## cRx(a,b,-0.2)
  #rx q1,-0.05
  #cz q0,q1
  #rx q1,0.1
  #cz q0,q1
  #rx q1,-0.05
  #display
  
  # Rx(b,-0.2)
  rx q1,-0.2
  display

  # init
  measure
  display