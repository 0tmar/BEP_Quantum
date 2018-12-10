#function: exp(2*pi*i*A/16) test
# 
# The matrix A is: (1/4)*[15,9,5,-3;9,15,3,-5;5,3,15,-9;-3,-5,-9,15];
# 
# The eigenvectors are: [-1;1;1;1], [1;-1;1;1], [1;1;-1;1], [1;1;1;-1] with eigenvalues 1, 2, 4, 8
# 
# The value of the first item in one of the vectors is represented by the amplitude of the |00> state, 
# The value of the second item by the amplitude of the |01> state, 
# The value of the third item by the amplitude of the |10> state, 
# The value of the fourth item by the amplitude of the |11> state.
# 
# The eigenvectors can be built with q1 and q2 using:
#    h q1
#    h q2
#    (x q1)
#    [x q2]
#    cz q1,q2
#    [x q2]
#    (x q1)
# The commands in brackets are optional:
# - if none of the bracketed commands are performed, we get [1;1;1;-1]
# - If only the optional commands in the square brackets are performed, we get [1;1;-1;1]
# - If only the optional commands in the round brackets are performed, we get [1;-1;1;1]
# - If both the optional commands in the round and quare brackets are performed, we get [-1;1;1;1]

qubits 3

.init
  x q0
  h q1
  h q2
  #x q1
  #x q2
  cz q1,q2
  #x q2
  #x q1
  display

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/16) on q1 and q2, controlled by q0
  
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
  
  # cRx(a,c,0.2)
  rx q2,0.05
  cz q0,q2
  rx q2,-0.1
  cz q0,q2
  rx q2,0.05
  
  # cVdag(a,c)
  t q0
  x q0
  tdag q0
  x q0
  h q2
  tdag q2
  tdag q2
  cx q0,q2
  h q2
  t q2
  cx q0,q2
  t q2
  h q2
  s q2
  x q2
  cz q0,q2
  t q0
  t q2
  cx q0,q2
  tdag q2
  cx q0,q2
  t q0
  x q0
  tdag q0
  x q0
  h q2
  tdag q2
  tdag q2
  cx q0,q2
  h q2
  t q2
  cx q0,q2
  t q2
  h q2
  s q2
  x q2
  
  # Rz(a,0.38)
  rz q0,0.38
  
  # cRx(a,b,0.98)
  rx q1,0.245
  cz q0,q1
  rx q1,-0.49
  cz q0,q1
  rx q1,0.245
  
  # Rz(a,1.88)
  rz q0,1.88
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
  # cRx(a,b,0.59)
  rx q1,0.1475
  cz q0,q1
  rx q1,-0.295
  cz q0,q1
  rx q1,0.1475
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
  # cZ(a,c)
  cz q0,q2
  
  display

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/16) on q1 and q2, controlled by q0
  
  # cZ(a,c)
  cz q0,q2
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
  # cRx(a,b,-0.59)
  rx q1,-0.1475
  cz q0,q1
  rx q1,0.295
  cz q0,q1
  rx q1,-0.1475
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
  # Rz(a,-1.88)
  rz q0,-1.88
  
  # cRx(a,b,-0.98)
  rx q1,-0.245
  cz q0,q1
  rx q1,0.49
  cz q0,q1
  rx q1,-0.245
  
  # Rz(a,-0.38)
  rz q0,-0.38
  
  # cV(a,c)
  t q0
  x q0
  tdag q0
  x q0
  h q2
  tdag q2
  tdag q2
  cx q0,q2
  h q2
  t q2
  cx q0,q2
  t q2
  h q2
  s q2
  x q2
  t q0
  t q2
  cx q0,q2
  tdag q2
  cx q0,q2
  t q0
  x q0
  tdag q0
  x q0
  h q2
  tdag q2
  tdag q2
  cx q0,q2
  h q2
  t q2
  cx q0,q2
  t q2
  h q2
  s q2
  x q2
  
  # cRx(a,c,-0.2)
  rx q2,-0.05
  cz q0,q2
  rx q2,0.1
  cz q0,q2
  rx q2,-0.05
  
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

.result
  measure
  display