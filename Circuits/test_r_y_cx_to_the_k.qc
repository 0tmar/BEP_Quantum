#function: Ry(c*(x^k)) rotation for c=1.1 and k=7

qubits 11

.init
  map q0,c0
  map q1,c1
  map q2,c2
  map q3,c3
  map q4,c4
  map q5,a0
  map q6,a1
  map q7,a2
  map q8,a3
  map q9,a4
  map q10,r
  x c1
  x c2
  x c3
  x c4
  display

.ry_c_x_to_the_k
  ##### Performs Ry(c*(x^k)) for c=1.1 and k=7
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 1.1 and theta = 1.0, controlled by [0]
  ry r,0.55
  cx c0,r
  ry r,-1.1
  cx c0,r
  ry r,0.55
  
  # Rotation c*theta with c = 1.1 and theta = 0.0078125, controlled by [1]
  ry r,0.004296875
  cx c1,r
  ry r,-0.00859375
  cx c1,r
  ry r,0.004296875
  
  # Rotation c*theta with c = 1.1 and theta = 6.103515625e-05, controlled by [2]
  ry r,3.35693359375e-05
  cx c2,r
  ry r,-6.7138671875e-05
  cx c2,r
  ry r,3.35693359375e-05
  
  # Rotation c*theta with c = 1.1 and theta = 4.76837158203125e-07, controlled by [3]
  ry r,2.6226043701171877e-07
  cx c3,r
  ry r,-5.245208740234375e-07
  cx c3,r
  ry r,2.6226043701171877e-07
  
  # Rotation c*theta with c = 1.1 and theta = 3.725290298461914e-09, controlled by [4]
  ry r,2.048909664154053e-09
  cx c4,r
  ry r,-4.097819328308106e-09
  cx c4,r
  ry r,2.048909664154053e-09
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 1.1 and theta = 16.078125, controlled by [0, 1]
  ry r,8.84296875
  cx c1,r
  ry r,-8.84296875
  cx c0,c1
  cx c0,r
  ry r,8.84296875
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-8.84296875
  cx c0,r
  
  # Rotation c*theta with c = 1.1 and theta = 0.1256103515625, controlled by [1, 2]
  ry r,0.069085693359375
  cx c2,r
  ry r,-0.069085693359375
  cx c1,c2
  cx c1,r
  ry r,0.069085693359375
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.069085693359375
  cx c1,r
  
  # Rotation c*theta with c = 1.1 and theta = 0.029439926147460938, controlled by [1, 3]
  ry r,0.016191959381103516
  cx c3,r
  ry r,-0.016191959381103516
  cx c1,c3
  cx c1,r
  ry r,0.016191959381103516
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.016191959381103516
  cx c1,r
  
  # Rotation c*theta with c = 1.1 and theta = 0.010005444288253784, controlled by [1, 4]
  ry r,0.0055029943585395815
  cx c4,r
  ry r,-0.0055029943585395815
  cx c1,c4
  cx c1,r
  ry r,0.0055029943585395815
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.0055029943585395815
  cx c1,r
  
  # Rotation c*theta with c = 1.1 and theta = 0.00022999942302703857, controlled by [2, 4]
  ry r,0.00012649968266487122
  cx c4,r
  ry r,-0.00012649968266487122
  cx c2,c4
  cx c2,r
  ry r,0.00012649968266487122
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-0.00012649968266487122
  cx c2,r
  
  # Rotation c*theta with c = 1.1 and theta = 0.0009813308715820312, controlled by [2, 3]
  ry r,0.0005397319793701172
  cx c3,r
  ry r,-0.0005397319793701172
  cx c2,c3
  cx c2,r
  ry r,0.0005397319793701172
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-0.0005397319793701172
  cx c2,r
  
  # Rotation c*theta with c = 1.1 and theta = 0.5286306738853455, controlled by [0, 4]
  ry r,0.29074687063694005
  cx c4,r
  ry r,-0.29074687063694005
  cx c0,c4
  cx c0,r
  ry r,0.29074687063694005
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.29074687063694005
  cx c0,r
  
  # Rotation c*theta with c = 1.1 and theta = 1.2806968688964844, controlled by [0, 3]
  ry r,0.7043832778930664
  cx c3,r
  ry r,-0.7043832778930664
  cx c0,c3
  cx c0,r
  ry r,0.7043832778930664
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.7043832778930664
  cx c0,r
  
  # Rotation c*theta with c = 1.1 and theta = 7.666647434234619e-06, controlled by [3, 4]
  ry r,4.216656088829041e-06
  cx c4,r
  ry r,-4.216656088829041e-06
  cx c3,c4
  cx c3,r
  ry r,4.216656088829041e-06
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-4.216656088829041e-06
  cx c3,r
  
  # Rotation c*theta with c = 1.1 and theta = 3.768310546875, controlled by [0, 2]
  ry r,2.07257080078125
  cx c2,r
  ry r,-2.07257080078125
  cx c0,c2
  cx c0,r
  ry r,2.07257080078125
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-2.07257080078125
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 1.1 and theta = 0.5205960273742676, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.2863278150558472
  cx a0,r
  ry r,-0.2863278150558472
  cx c4,a0
  cx c4,r
  ry r,0.2863278150558472
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.2863278150558472
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 1.1 and theta = 29.28515625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,16.1068359375
  cx a0,r
  ry r,-16.1068359375
  cx c2,a0
  cx c2,r
  ry r,16.1068359375
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-16.1068359375
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 1.1 and theta = 3.242156982421875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,1.7831863403320314
  cx a0,r
  ry r,-1.7831863403320314
  cx c3,a0
  cx c3,r
  ry r,1.7831863403320314
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-1.7831863403320314
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 1.1 and theta = 11.524749755859375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,6.338612365722657
  cx a0,r
  ry r,-6.338612365722657
  cx c3,a0
  cx c3,r
  ry r,6.338612365722657
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-6.338612365722657
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 1.1 and theta = 0.228790283203125, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.12583465576171876
  cx a0,r
  ry r,-0.12583465576171876
  cx c3,a0
  cx c3,r
  ry r,0.12583465576171876
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.12583465576171876
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 1.1 and theta = 0.0253293514251709, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.013931143283843995
  cx a0,r
  ry r,-0.013931143283843995
  cx c4,a0
  cx c4,r
  ry r,0.013931143283843995
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.013931143283843995
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 1.1 and theta = 5.112793922424316, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,2.8120366573333744
  cx a0,r
  ry r,-2.8120366573333744
  cx c4,a0
  cx c4,r
  ry r,2.8120366573333744
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-2.8120366573333744
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 1.1 and theta = 0.09003710746765137, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.04952040910720826
  cx a0,r
  ry r,-0.04952040910720826
  cx c4,a0
  cx c4,r
  ry r,0.04952040910720826
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.04952040910720826
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 1.1 and theta = 0.001787424087524414, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,0.0009830832481384278
  cx a0,r
  ry r,-0.0009830832481384278
  cx c4,a0
  cx c4,r
  ry r,0.0009830832481384278
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0009830832481384278
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 1.1 and theta = 1.4123454093933105, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.7767899751663209
  cx a0,r
  ry r,-0.7767899751663209
  cx c4,a0
  cx c4,r
  ry r,0.7767899751663209
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.7767899751663209
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 1.1 and theta = 0.11640787124633789, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.06402432918548584
  cx a1,r
  ry r,-0.06402432918548584
  cx c4,a1
  cx c4,r
  ry r,0.06402432918548584
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.06402432918548584
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 1.1 and theta = 2.8496646881103516, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,1.5673155784606936
  cx a1,r
  ry r,-1.5673155784606936
  cx c4,a1
  cx c4,r
  ry r,1.5673155784606936
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-1.5673155784606936
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 1.1 and theta = 6.841678619384766, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,3.7629232406616215
  cx a1,r
  ry r,-3.7629232406616215
  cx c4,a1
  cx c4,r
  ry r,3.7629232406616215
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-3.7629232406616215
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 1.1 and theta = 0.9281587600708008, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,0.5104873180389404
  cx a1,r
  ry r,-0.5104873180389404
  cx c4,a1
  cx c4,r
  ry r,0.5104873180389404
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.5104873180389404
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 1.1 and theta = 14.90020751953125, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,8.195114135742188
  cx a1,r
  ry r,-8.195114135742188
  cx c3,a1
  cx c3,r
  ry r,8.195114135742188
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-8.195114135742188
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 1.1 and theta = 2.58270263671875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,1.4204864501953127
  cx a2,r
  ry r,-1.4204864501953127
  cx c4,a2
  cx c4,r
  ry r,1.4204864501953127
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-1.4204864501953127
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.result
  display