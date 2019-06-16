#function: Ry(c*(x^k)) rotation for c=10.0 and k=7

qubits 9

.init
  map q0,c0
  map q1,c1
  map q2,c2
  map q3,c3
  map q4,c4
  map q5,a0
  map q6,a1
  map q7,a2
  map q8,r
  x c1
  x c3
  x c4
  display

.ry_c_x_to_the_k
  ##### Performs Ry(c*(x^k)) for c=10.0 and k=7
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 10.0 and theta = 1.0, controlled by [0]
  ry r,5.0
  cx c0,r
  ry r,-10.0
  cx c0,r
  ry r,5.0
  
  # Rotation c*theta with c = 10.0 and theta = 0.0078125, controlled by [1]
  ry r,0.0390625
  cx c1,r
  ry r,-0.078125
  cx c1,r
  ry r,0.0390625
  
  # Rotation c*theta with c = 10.0 and theta = 6.103515625e-05, controlled by [2]
  ry r,0.00030517578125
  cx c2,r
  ry r,-0.0006103515625
  cx c2,r
  ry r,0.00030517578125
  
  # Rotation c*theta with c = 10.0 and theta = 4.76837158203125e-07, controlled by [3]
  ry r,2.384185791015625e-06
  cx c3,r
  ry r,-4.76837158203125e-06
  cx c3,r
  ry r,2.384185791015625e-06
  
  # Rotation c*theta with c = 10.0 and theta = 3.725290298461914e-09, controlled by [4]
  ry r,1.862645149230957e-08
  cx c4,r
  ry r,-3.725290298461914e-08
  cx c4,r
  ry r,1.862645149230957e-08
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 10.0 and theta = 16.078125, controlled by [0, 1]
  ry r,80.390625
  cx c1,r
  ry r,-80.390625
  cx c0,c1
  cx c0,r
  ry r,80.390625
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-80.390625
  cx c0,r
  
  # Rotation c*theta with c = 10.0 and theta = 0.1256103515625, controlled by [1, 2]
  ry r,0.6280517578125
  cx c2,r
  ry r,-0.6280517578125
  cx c1,c2
  cx c1,r
  ry r,0.6280517578125
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.6280517578125
  cx c1,r
  
  # Rotation c*theta with c = 10.0 and theta = 0.029439926147460938, controlled by [1, 3]
  ry r,0.1471996307373047
  cx c3,r
  ry r,-0.1471996307373047
  cx c1,c3
  cx c1,r
  ry r,0.1471996307373047
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.1471996307373047
  cx c1,r
  
  # Rotation c*theta with c = 10.0 and theta = 0.010005444288253784, controlled by [1, 4]
  ry r,0.05002722144126892
  cx c4,r
  ry r,-0.05002722144126892
  cx c1,c4
  cx c1,r
  ry r,0.05002722144126892
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.05002722144126892
  cx c1,r
  
  # Rotation c*theta with c = 10.0 and theta = 0.00022999942302703857, controlled by [2, 4]
  ry r,0.0011499971151351929
  cx c4,r
  ry r,-0.0011499971151351929
  cx c2,c4
  cx c2,r
  ry r,0.0011499971151351929
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-0.0011499971151351929
  cx c2,r
  
  # Rotation c*theta with c = 10.0 and theta = 0.0009813308715820312, controlled by [2, 3]
  ry r,0.004906654357910156
  cx c3,r
  ry r,-0.004906654357910156
  cx c2,c3
  cx c2,r
  ry r,0.004906654357910156
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-0.004906654357910156
  cx c2,r
  
  # Rotation c*theta with c = 10.0 and theta = 0.5286306738853455, controlled by [0, 4]
  ry r,2.6431533694267273
  cx c4,r
  ry r,-2.6431533694267273
  cx c0,c4
  cx c0,r
  ry r,2.6431533694267273
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-2.6431533694267273
  cx c0,r
  
  # Rotation c*theta with c = 10.0 and theta = 1.2806968688964844, controlled by [0, 3]
  ry r,6.403484344482422
  cx c3,r
  ry r,-6.403484344482422
  cx c0,c3
  cx c0,r
  ry r,6.403484344482422
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-6.403484344482422
  cx c0,r
  
  # Rotation c*theta with c = 10.0 and theta = 7.666647434234619e-06, controlled by [3, 4]
  ry r,3.8333237171173096e-05
  cx c4,r
  ry r,-3.8333237171173096e-05
  cx c3,c4
  cx c3,r
  ry r,3.8333237171173096e-05
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-3.8333237171173096e-05
  cx c3,r
  
  # Rotation c*theta with c = 10.0 and theta = 3.768310546875, controlled by [0, 2]
  ry r,18.841552734375
  cx c2,r
  ry r,-18.841552734375
  cx c0,c2
  cx c0,r
  ry r,18.841552734375
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-18.841552734375
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 10.0 and theta = 0.5205960273742676, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,2.602980136871338
  cx a0,r
  ry r,-2.602980136871338
  cx c4,a0
  cx c4,r
  ry r,2.602980136871338
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-2.602980136871338
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 10.0 and theta = 29.28515625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,146.42578125
  cx a0,r
  ry r,-146.42578125
  cx c2,a0
  cx c2,r
  ry r,146.42578125
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-146.42578125
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 10.0 and theta = 3.242156982421875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,16.210784912109375
  cx a0,r
  ry r,-16.210784912109375
  cx c3,a0
  cx c3,r
  ry r,16.210784912109375
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-16.210784912109375
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 10.0 and theta = 11.524749755859375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,57.623748779296875
  cx a0,r
  ry r,-57.623748779296875
  cx c3,a0
  cx c3,r
  ry r,57.623748779296875
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-57.623748779296875
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 10.0 and theta = 0.228790283203125, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,1.143951416015625
  cx a0,r
  ry r,-1.143951416015625
  cx c3,a0
  cx c3,r
  ry r,1.143951416015625
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-1.143951416015625
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 10.0 and theta = 0.0253293514251709, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.1266467571258545
  cx a0,r
  ry r,-0.1266467571258545
  cx c4,a0
  cx c4,r
  ry r,0.1266467571258545
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.1266467571258545
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 10.0 and theta = 5.112793922424316, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,25.563969612121582
  cx a0,r
  ry r,-25.563969612121582
  cx c4,a0
  cx c4,r
  ry r,25.563969612121582
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-25.563969612121582
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 10.0 and theta = 0.09003710746765137, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.45018553733825684
  cx a0,r
  ry r,-0.45018553733825684
  cx c4,a0
  cx c4,r
  ry r,0.45018553733825684
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.45018553733825684
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 10.0 and theta = 0.001787424087524414, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,0.00893712043762207
  cx a0,r
  ry r,-0.00893712043762207
  cx c4,a0
  cx c4,r
  ry r,0.00893712043762207
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.00893712043762207
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 10.0 and theta = 1.4123454093933105, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,7.061727046966553
  cx a0,r
  ry r,-7.061727046966553
  cx c4,a0
  cx c4,r
  ry r,7.061727046966553
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-7.061727046966553
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 10.0 and theta = 0.11640787124633789, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.5820393562316895
  cx a1,r
  ry r,-0.5820393562316895
  cx c4,a1
  cx c4,r
  ry r,0.5820393562316895
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.5820393562316895
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 10.0 and theta = 2.8496646881103516, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,14.248323440551758
  cx a1,r
  ry r,-14.248323440551758
  cx c4,a1
  cx c4,r
  ry r,14.248323440551758
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-14.248323440551758
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 10.0 and theta = 6.841678619384766, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,34.20839309692383
  cx a1,r
  ry r,-34.20839309692383
  cx c4,a1
  cx c4,r
  ry r,34.20839309692383
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-34.20839309692383
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 10.0 and theta = 0.9281587600708008, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,4.640793800354004
  cx a1,r
  ry r,-4.640793800354004
  cx c4,a1
  cx c4,r
  ry r,4.640793800354004
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-4.640793800354004
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 10.0 and theta = 14.90020751953125, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,74.50103759765625
  cx a1,r
  ry r,-74.50103759765625
  cx c3,a1
  cx c3,r
  ry r,74.50103759765625
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-74.50103759765625
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 10.0 and theta = 2.58270263671875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,12.91351318359375
  cx a2,r
  ry r,-12.91351318359375
  cx c4,a2
  cx c4,r
  ry r,12.91351318359375
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-12.91351318359375
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.result
  display