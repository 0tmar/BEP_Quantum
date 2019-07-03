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
# - if none of the bracketed commands are performed, we get [1;1;1;-1], and a rotation of pi/8
# - If only the optional commands in the round brackets are performed, we get [1;-1;1;1], and a rotation of pi/4
# - If only the optional commands in the square brackets are performed, we get [1;1;-1;1], and a rotation of pi/2
# - If both the optional commands in the round and quare brackets are performed, we get [-1;1;1;1], and a rotation of pi

qubits 4

.init
  x q0
  h q1
  h q2
  #x q1
  #x q2
  cz q1,q2
  #x q2
  #x q1
  x q3
  display

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(2*pi*i*A/2) on q1 and q2, controlled by q0
  
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
  
  # cRx(a,c,-1.568)
  rx q2,-0.392
  cz q0,q2
  rx q2,0.784
  cz q0,q2
  rx q2,-0.392
  
  # Rz(a,3.0)
  rz q0,3.0
  
  # cRx(a,b,-7.856)
  rx q1,-1.964
  cz q0,q1
  rx q1,3.928
  cz q0,q1
  rx q1,-1.964
  
  # Rz(a,15.064)
  rz q0,15.064
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
  # cRx(a,b,-4.712)
  rx q1,-1.178
  cz q0,q1
  rx q1,2.356
  cz q0,q1
  rx q1,-1.178
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
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
  
  # Rz(d,18.064)
  rz q3,18.064
  
  display

.expA
  # performs exp(sign*i*A*t0*2^(n)/16) = exp(-2*pi*i*A/2) on q1 and q2, controlled by q0
  
  # Rz(d,-18.064)
  rz q3,-18.064
  
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
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
  # cRx(a,b,4.712)
  rx q1,1.178
  cz q0,q1
  rx q1,-2.356
  cz q0,q1
  rx q1,1.178
  
  # Toffoli(a,b,c)
  toffoli q0,q1,q2
  
  # Rz(a,-15.064)
  rz q0,-15.064
  
  # cRx(a,b,7.856)
  rx q1,1.964
  cz q0,q1
  rx q1,-3.928
  cz q0,q1
  rx q1,1.964
  
  # Rz(a,-3.0)
  rz q0,-3.0
  
  # cRx(a,c,1.568)
  rx q2,0.392
  cz q0,q2
  rx q2,-0.784
  cz q0,q2
  rx q2,0.392
  
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

.result
  display
  measure
  display