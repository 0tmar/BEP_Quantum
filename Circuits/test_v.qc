#function: V-Vdag test

qubits 1

.test

  # init
  x q0
  display
  
  # V
  h q0
  s q0
  h q0
  display
  
  # Vdag
  h q0
  z q0
  s q0
  h q0
  display

  # meas
  measure
  display