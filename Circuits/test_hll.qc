#function: Eigenvalue inversion test

qubits 10

.init
  map q0,l0
  map q1,l1
  map q2,m0
  map q3,m1
  map q4,m2
  map q5,m3
  map q6,m4
  map q7,m5
  map q8,c0
  map q9,c1
  h l0
  h l1
  h m0
  h m1
  h m2
  h m3
  h m4
  h m5
  x c1
  display

.cRzz
  # Performs cRzz(2*pi/2^(m-i)) on L, controlled by the i-th qubit of M, for all i in M
  rz m0,0.09817477042468103
  rz m1,0.19634954084936207
  rz m2,0.39269908169872414
  rz m3,0.7853981633974483
  rz m4,1.5707963267948966
  rz m5,3.141592653589793

.exp_iH0t0
  
  ####### Performs exp(i*H0*t0) with t0 = 2*pi
  
  ##### Performs G(l - kl) with l = 2 and kl = 0
  
  ### Performs exp(i*H0*t/2^m) with t = 25.132741228718345
  
  # ccRz(l0,c0,m0,12.566370614359172)
  rz m0,3.141592653589793
  cx c0,m0
  rz m0,-3.141592653589793
  cx l0,c0
  cx l0,m0
  rz m0,3.141592653589793
  cx c0,m0
  cx l0,c0
  cx l0,m0
  rz m0,-3.141592653589793
  cx l0,m0
  
  # ccRz(l0,c0,m1,6.283185307179586)
  rz m1,1.5707963267948966
  cx c0,m1
  rz m1,-1.5707963267948966
  cx l0,c0
  cx l0,m1
  rz m1,1.5707963267948966
  cx c0,m1
  cx l0,c0
  cx l0,m1
  rz m1,-1.5707963267948966
  cx l0,m1
  
  # ccRz(l0,c0,m2,3.141592653589793)
  rz m2,0.7853981633974483
  cx c0,m2
  rz m2,-0.7853981633974483
  cx l0,c0
  cx l0,m2
  rz m2,0.7853981633974483
  cx c0,m2
  cx l0,c0
  cx l0,m2
  rz m2,-0.7853981633974483
  cx l0,m2
  
  # ccRz(l0,c0,m3,1.5707963267948966)
  rz m3,0.39269908169872414
  cx c0,m3
  rz m3,-0.39269908169872414
  cx l0,c0
  cx l0,m3
  rz m3,0.39269908169872414
  cx c0,m3
  cx l0,c0
  cx l0,m3
  rz m3,-0.39269908169872414
  cx l0,m3
  
  # ccRz(l0,c0,m4,0.7853981633974483)
  rz m4,0.19634954084936207
  cx c0,m4
  rz m4,-0.19634954084936207
  cx l0,c0
  cx l0,m4
  rz m4,0.19634954084936207
  cx c0,m4
  cx l0,c0
  cx l0,m4
  rz m4,-0.19634954084936207
  cx l0,m4
  
  # ccRz(l0,c0,m5,0.39269908169872414)
  rz m5,0.09817477042468103
  cx c0,m5
  rz m5,-0.09817477042468103
  cx l0,c0
  cx l0,m5
  rz m5,0.09817477042468103
  cx c0,m5
  cx l0,c0
  cx l0,m5
  rz m5,-0.09817477042468103
  cx l0,m5
  
  ### Performs exp(i*H0*t/2^m) with t = 50.26548245743669
  
  # ccRz(l0,c1,m0,25.132741228718345)
  rz m0,6.283185307179586
  cx c1,m0
  rz m0,-6.283185307179586
  cx l0,c1
  cx l0,m0
  rz m0,6.283185307179586
  cx c1,m0
  cx l0,c1
  cx l0,m0
  rz m0,-6.283185307179586
  cx l0,m0
  
  # ccRz(l0,c1,m1,12.566370614359172)
  rz m1,3.141592653589793
  cx c1,m1
  rz m1,-3.141592653589793
  cx l0,c1
  cx l0,m1
  rz m1,3.141592653589793
  cx c1,m1
  cx l0,c1
  cx l0,m1
  rz m1,-3.141592653589793
  cx l0,m1
  
  # ccRz(l0,c1,m2,6.283185307179586)
  rz m2,1.5707963267948966
  cx c1,m2
  rz m2,-1.5707963267948966
  cx l0,c1
  cx l0,m2
  rz m2,1.5707963267948966
  cx c1,m2
  cx l0,c1
  cx l0,m2
  rz m2,-1.5707963267948966
  cx l0,m2
  
  # ccRz(l0,c1,m3,3.141592653589793)
  rz m3,0.7853981633974483
  cx c1,m3
  rz m3,-0.7853981633974483
  cx l0,c1
  cx l0,m3
  rz m3,0.7853981633974483
  cx c1,m3
  cx l0,c1
  cx l0,m3
  rz m3,-0.7853981633974483
  cx l0,m3
  
  # ccRz(l0,c1,m4,1.5707963267948966)
  rz m4,0.39269908169872414
  cx c1,m4
  rz m4,-0.39269908169872414
  cx l0,c1
  cx l0,m4
  rz m4,0.39269908169872414
  cx c1,m4
  cx l0,c1
  cx l0,m4
  rz m4,-0.39269908169872414
  cx l0,m4
  
  # ccRz(l0,c1,m5,0.7853981633974483)
  rz m5,0.19634954084936207
  cx c1,m5
  rz m5,-0.19634954084936207
  cx l0,c1
  cx l0,m5
  rz m5,0.19634954084936207
  cx c1,m5
  cx l0,c1
  cx l0,m5
  rz m5,-0.19634954084936207
  cx l0,m5
  
  ##### Performs G(l - kl) with l = 2 and kl = 1
  
  ### Performs exp(i*H0*t/2^m) with t = 50.26548245743669
  
  # ccRz(l1,c0,m0,25.132741228718345)
  rz m0,6.283185307179586
  cx c0,m0
  rz m0,-6.283185307179586
  cx l1,c0
  cx l1,m0
  rz m0,6.283185307179586
  cx c0,m0
  cx l1,c0
  cx l1,m0
  rz m0,-6.283185307179586
  cx l1,m0
  
  # ccRz(l1,c0,m1,12.566370614359172)
  rz m1,3.141592653589793
  cx c0,m1
  rz m1,-3.141592653589793
  cx l1,c0
  cx l1,m1
  rz m1,3.141592653589793
  cx c0,m1
  cx l1,c0
  cx l1,m1
  rz m1,-3.141592653589793
  cx l1,m1
  
  # ccRz(l1,c0,m2,6.283185307179586)
  rz m2,1.5707963267948966
  cx c0,m2
  rz m2,-1.5707963267948966
  cx l1,c0
  cx l1,m2
  rz m2,1.5707963267948966
  cx c0,m2
  cx l1,c0
  cx l1,m2
  rz m2,-1.5707963267948966
  cx l1,m2
  
  # ccRz(l1,c0,m3,3.141592653589793)
  rz m3,0.7853981633974483
  cx c0,m3
  rz m3,-0.7853981633974483
  cx l1,c0
  cx l1,m3
  rz m3,0.7853981633974483
  cx c0,m3
  cx l1,c0
  cx l1,m3
  rz m3,-0.7853981633974483
  cx l1,m3
  
  # ccRz(l1,c0,m4,1.5707963267948966)
  rz m4,0.39269908169872414
  cx c0,m4
  rz m4,-0.39269908169872414
  cx l1,c0
  cx l1,m4
  rz m4,0.39269908169872414
  cx c0,m4
  cx l1,c0
  cx l1,m4
  rz m4,-0.39269908169872414
  cx l1,m4
  
  # ccRz(l1,c0,m5,0.7853981633974483)
  rz m5,0.19634954084936207
  cx c0,m5
  rz m5,-0.19634954084936207
  cx l1,c0
  cx l1,m5
  rz m5,0.19634954084936207
  cx c0,m5
  cx l1,c0
  cx l1,m5
  rz m5,-0.19634954084936207
  cx l1,m5
  
  ### Performs exp(i*H0*t/2^m) with t = 100.53096491487338
  
  # ccRz(l1,c1,m0,50.26548245743669)
  rz m0,12.566370614359172
  cx c1,m0
  rz m0,-12.566370614359172
  cx l1,c1
  cx l1,m0
  rz m0,12.566370614359172
  cx c1,m0
  cx l1,c1
  cx l1,m0
  rz m0,-12.566370614359172
  cx l1,m0
  
  # ccRz(l1,c1,m1,25.132741228718345)
  rz m1,6.283185307179586
  cx c1,m1
  rz m1,-6.283185307179586
  cx l1,c1
  cx l1,m1
  rz m1,6.283185307179586
  cx c1,m1
  cx l1,c1
  cx l1,m1
  rz m1,-6.283185307179586
  cx l1,m1
  
  # ccRz(l1,c1,m2,12.566370614359172)
  rz m2,3.141592653589793
  cx c1,m2
  rz m2,-3.141592653589793
  cx l1,c1
  cx l1,m2
  rz m2,3.141592653589793
  cx c1,m2
  cx l1,c1
  cx l1,m2
  rz m2,-3.141592653589793
  cx l1,m2
  
  # ccRz(l1,c1,m3,6.283185307179586)
  rz m3,1.5707963267948966
  cx c1,m3
  rz m3,-1.5707963267948966
  cx l1,c1
  cx l1,m3
  rz m3,1.5707963267948966
  cx c1,m3
  cx l1,c1
  cx l1,m3
  rz m3,-1.5707963267948966
  cx l1,m3
  
  # ccRz(l1,c1,m4,3.141592653589793)
  rz m4,0.7853981633974483
  cx c1,m4
  rz m4,-0.7853981633974483
  cx l1,c1
  cx l1,m4
  rz m4,0.7853981633974483
  cx c1,m4
  cx l1,c1
  cx l1,m4
  rz m4,-0.7853981633974483
  cx l1,m4
  
  # ccRz(l1,c1,m5,1.5707963267948966)
  rz m5,0.39269908169872414
  cx c1,m5
  rz m5,-0.39269908169872414
  cx l1,c1
  cx l1,m5
  rz m5,0.39269908169872414
  cx c1,m5
  cx l1,c1
  cx l1,m5
  rz m5,-0.39269908169872414
  cx l1,m5

.measure
  display
  measure