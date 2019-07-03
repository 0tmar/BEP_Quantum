#function: Ancilla rotation for c=1.1, approximated using the m=4 first terms of the taylor expansion of arcsin(x)

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

.ry_c_x_to_the_p
  ##### Performs Ry(c*(x^p)) for c=1.1 and p=1
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 1.1 and theta = 1.0, controlled by [0]
  ry r,0.55
  cx c0,r
  ry r,-1.1
  cx c0,r
  ry r,0.55
  
  # Rotation c*theta with c = 1.1 and theta = 0.5, controlled by [1]
  ry r,0.275
  cx c1,r
  ry r,-0.55
  cx c1,r
  ry r,0.275
  
  # Rotation c*theta with c = 1.1 and theta = 0.25, controlled by [2]
  ry r,0.1375
  cx c2,r
  ry r,-0.275
  cx c2,r
  ry r,0.1375
  
  # Rotation c*theta with c = 1.1 and theta = 0.125, controlled by [3]
  ry r,0.06875
  cx c3,r
  ry r,-0.1375
  cx c3,r
  ry r,0.06875
  
  # Rotation c*theta with c = 1.1 and theta = 0.0625, controlled by [4]
  ry r,0.034375
  cx c4,r
  ry r,-0.06875
  cx c4,r
  ry r,0.034375
  

.ry_c_x_to_the_p
  ##### Performs Ry(c*(x^p)) for c=0.22183333333333338 and p=3
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 1.0, controlled by [0]
  ry r,0.11091666666666669
  cx c0,r
  ry r,-0.22183333333333338
  cx c0,r
  ry r,0.11091666666666669
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.125, controlled by [1]
  ry r,0.013864583333333336
  cx c1,r
  ry r,-0.027729166666666673
  cx c1,r
  ry r,0.013864583333333336
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.015625, controlled by [2]
  ry r,0.001733072916666667
  cx c2,r
  ry r,-0.003466145833333334
  cx c2,r
  ry r,0.001733072916666667
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.001953125, controlled by [3]
  ry r,0.00021663411458333338
  cx c3,r
  ry r,-0.00043326822916666676
  cx c3,r
  ry r,0.00021663411458333338
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.000244140625, controlled by [4]
  ry r,2.7079264322916673e-05
  cx c4,r
  ry r,-5.4158528645833345e-05
  cx c4,r
  ry r,2.7079264322916673e-05
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 2.25, controlled by [0, 1]
  ry r,0.24956250000000005
  cx c1,r
  ry r,-0.24956250000000005
  cx c0,c1
  cx c0,r
  ry r,0.24956250000000005
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.24956250000000005
  cx c0,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.28125, controlled by [1, 2]
  ry r,0.031195312500000006
  cx c2,r
  ry r,-0.031195312500000006
  cx c1,c2
  cx c1,r
  ry r,0.031195312500000006
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.031195312500000006
  cx c1,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.1171875, controlled by [1, 3]
  ry r,0.012998046875000002
  cx c3,r
  ry r,-0.012998046875000002
  cx c1,c3
  cx c1,r
  ry r,0.012998046875000002
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.012998046875000002
  cx c1,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.052734375, controlled by [1, 4]
  ry r,0.005849121093750001
  cx c4,r
  ry r,-0.005849121093750001
  cx c1,c4
  cx c1,r
  ry r,0.005849121093750001
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.005849121093750001
  cx c1,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.0146484375, controlled by [2, 4]
  ry r,0.0016247558593750003
  cx c4,r
  ry r,-0.0016247558593750003
  cx c2,c4
  cx c2,r
  ry r,0.0016247558593750003
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-0.0016247558593750003
  cx c2,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.03515625, controlled by [2, 3]
  ry r,0.0038994140625000008
  cx c3,r
  ry r,-0.0038994140625000008
  cx c2,c3
  cx c2,r
  ry r,0.0038994140625000008
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-0.0038994140625000008
  cx c2,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.19921875, controlled by [0, 4]
  ry r,0.022096679687500006
  cx c4,r
  ry r,-0.022096679687500006
  cx c0,c4
  cx c0,r
  ry r,0.022096679687500006
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.022096679687500006
  cx c0,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.421875, controlled by [0, 3]
  ry r,0.04679296875000001
  cx c3,r
  ry r,-0.04679296875000001
  cx c0,c3
  cx c0,r
  ry r,0.04679296875000001
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.04679296875000001
  cx c0,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.00439453125, controlled by [3, 4]
  ry r,0.0004874267578125001
  cx c4,r
  ry r,-0.0004874267578125001
  cx c3,c4
  cx c3,r
  ry r,0.0004874267578125001
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-0.0004874267578125001
  cx c3,r
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.9375, controlled by [0, 2]
  ry r,0.10398437500000002
  cx c2,r
  ry r,-0.10398437500000002
  cx c0,c2
  cx c0,r
  ry r,0.10398437500000002
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.10398437500000002
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.046875, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.005199218750000002
  cx a0,r
  ry r,-0.005199218750000002
  cx c4,a0
  cx c4,r
  ry r,0.005199218750000002
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.005199218750000002
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.75, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,0.08318750000000003
  cx a0,r
  ry r,-0.08318750000000003
  cx c2,a0
  cx c2,r
  ry r,0.08318750000000003
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-0.08318750000000003
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.1875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.020796875000000006
  cx a0,r
  ry r,-0.020796875000000006
  cx c3,a0
  cx c3,r
  ry r,0.020796875000000006
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.020796875000000006
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.04159375000000001
  cx a0,r
  ry r,-0.04159375000000001
  cx c3,a0
  cx c3,r
  ry r,0.04159375000000001
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.04159375000000001
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.09375, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.010398437500000003
  cx a0,r
  ry r,-0.010398437500000003
  cx c3,a0
  cx c3,r
  ry r,0.010398437500000003
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.010398437500000003
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.0234375, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.002599609375000001
  cx a0,r
  ry r,-0.002599609375000001
  cx c4,a0
  cx c4,r
  ry r,0.002599609375000001
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.002599609375000001
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.1875, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.020796875000000006
  cx a0,r
  ry r,-0.020796875000000006
  cx c4,a0
  cx c4,r
  ry r,0.020796875000000006
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.020796875000000006
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.046875, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.005199218750000002
  cx a0,r
  ry r,-0.005199218750000002
  cx c4,a0
  cx c4,r
  ry r,0.005199218750000002
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.005199218750000002
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.01171875, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,0.0012998046875000004
  cx a0,r
  ry r,-0.0012998046875000004
  cx c4,a0
  cx c4,r
  ry r,0.0012998046875000004
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0012998046875000004
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.22183333333333338 and theta = 0.09375, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.010398437500000003
  cx a0,r
  ry r,-0.010398437500000003
  cx c4,a0
  cx c4,r
  ry r,0.010398437500000003
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.010398437500000003
  cx c4,r
  toffoli c0,c2,a0
  

.ry_c_x_to_the_p
  ##### Performs Ry(c*(x^p)) for c=0.12078825000000004 and p=5
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 1.0, controlled by [0]
  ry r,0.06039412500000002
  cx c0,r
  ry r,-0.12078825000000004
  cx c0,r
  ry r,0.06039412500000002
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.03125, controlled by [1]
  ry r,0.0018873164062500006
  cx c1,r
  ry r,-0.0037746328125000013
  cx c1,r
  ry r,0.0018873164062500006
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.0009765625, controlled by [2]
  ry r,5.897863769531252e-05
  cx c2,r
  ry r,-0.00011795727539062504
  cx c2,r
  ry r,5.897863769531252e-05
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 3.0517578125e-05, controlled by [3]
  ry r,1.8430824279785163e-06
  cx c3,r
  ry r,-3.6861648559570325e-06
  cx c3,r
  ry r,1.8430824279785163e-06
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 9.5367431640625e-07, controlled by [4]
  ry r,5.759632587432863e-08
  cx c4,r
  ry r,-1.1519265174865727e-07
  cx c4,r
  ry r,5.759632587432863e-08
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 6.5625, controlled by [0, 1]
  ry r,0.3963364453125001
  cx c1,r
  ry r,-0.3963364453125001
  cx c0,c1
  cx c0,r
  ry r,0.3963364453125001
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.3963364453125001
  cx c0,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.205078125, controlled by [1, 2]
  ry r,0.012385513916015629
  cx c2,r
  ry r,-0.012385513916015629
  cx c1,c2
  cx c1,r
  ry r,0.012385513916015629
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.012385513916015629
  cx c1,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.0640869140625, controlled by [1, 3]
  ry r,0.003870473098754884
  cx c3,r
  ry r,-0.003870473098754884
  cx c1,c3
  cx c1,r
  ry r,0.003870473098754884
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.003870473098754884
  cx c1,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.02506256103515625, controlled by [1, 4]
  ry r,0.0015136314439773565
  cx c4,r
  ry r,-0.0015136314439773565
  cx c1,c4
  cx c1,r
  ry r,0.0015136314439773565
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.0015136314439773565
  cx c1,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.002002716064453125, controlled by [2, 4]
  ry r,0.00012095228433609013
  cx c4,r
  ry r,-0.00012095228433609013
  cx c2,c4
  cx c2,r
  ry r,0.00012095228433609013
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-0.00012095228433609013
  cx c2,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.00640869140625, controlled by [2, 3]
  ry r,0.0003870473098754884
  cx c3,r
  ry r,-0.0003870473098754884
  cx c2,c3
  cx c2,r
  ry r,0.0003870473098754884
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-0.0003870473098754884
  cx c2,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.3540802001953125, controlled by [0, 4]
  ry r,0.021384363870620735
  cx c4,r
  ry r,-0.021384363870620735
  cx c0,c4
  cx c0,r
  ry r,0.021384363870620735
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.021384363870620735
  cx c0,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.802001953125, controlled by [0, 3]
  ry r,0.04843620620727541
  cx c3,r
  ry r,-0.04843620620727541
  cx c0,c3
  cx c0,r
  ry r,0.04843620620727541
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.04843620620727541
  cx c0,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.0002002716064453125, controlled by [3, 4]
  ry r,1.2095228433609012e-05
  cx c4,r
  ry r,-1.2095228433609012e-05
  cx c3,c4
  cx c3,r
  ry r,1.2095228433609012e-05
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-1.2095228433609012e-05
  cx c3,r
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 2.05078125, controlled by [0, 2]
  ry r,0.1238551391601563
  cx c2,r
  ry r,-0.1238551391601563
  cx c0,c2
  cx c0,r
  ry r,0.1238551391601563
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.1238551391601563
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.205078125, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.012385513916015629
  cx a0,r
  ry r,-0.012385513916015629
  cx c4,a0
  cx c4,r
  ry r,0.012385513916015629
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.012385513916015629
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 6.5625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,0.3963364453125001
  cx a0,r
  ry r,-0.3963364453125001
  cx c2,a0
  cx c2,r
  ry r,0.3963364453125001
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-0.3963364453125001
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 1.0546875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.06369692871093752
  cx a0,r
  ry r,-0.06369692871093752
  cx c3,a0
  cx c3,r
  ry r,0.06369692871093752
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.06369692871093752
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 2.87109375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.1733971948242188
  cx a0,r
  ry r,-0.1733971948242188
  cx c3,a0
  cx c3,r
  ry r,0.1733971948242188
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.1733971948242188
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.205078125, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.012385513916015629
  cx a0,r
  ry r,-0.012385513916015629
  cx c3,a0
  cx c3,r
  ry r,0.012385513916015629
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.012385513916015629
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.032958984375, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.0019905290222167974
  cx a0,r
  ry r,-0.0019905290222167974
  cx c4,a0
  cx c4,r
  ry r,0.0019905290222167974
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0019905290222167974
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 1.34033203125, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.08094818023681644
  cx a0,r
  ry r,-0.08094818023681644
  cx c4,a0
  cx c4,r
  ry r,0.08094818023681644
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.08094818023681644
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.0897216796875, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.005418662338256838
  cx a0,r
  ry r,-0.005418662338256838
  cx c4,a0
  cx c4,r
  ry r,0.005418662338256838
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.005418662338256838
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.00640869140625, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,0.0003870473098754884
  cx a0,r
  ry r,-0.0003870473098754884
  cx c4,a0
  cx c4,r
  ry r,0.0003870473098754884
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0003870473098754884
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.487060546875, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.02941559555053712
  cx a0,r
  ry r,-0.02941559555053712
  cx c4,a0
  cx c4,r
  ry r,0.02941559555053712
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.02941559555053712
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.054931640625, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.0033175483703613293
  cx a1,r
  ry r,-0.0033175483703613293
  cx c4,a1
  cx c4,r
  ry r,0.0033175483703613293
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.0033175483703613293
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.3955078125, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,0.023886348266601572
  cx a1,r
  ry r,-0.023886348266601572
  cx c4,a1
  cx c4,r
  ry r,0.023886348266601572
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.023886348266601572
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.849609375, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.05131141479492189
  cx a1,r
  ry r,-0.05131141479492189
  cx c4,a1
  cx c4,r
  ry r,0.05131141479492189
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.05131141479492189
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.16845703125, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,0.01017381500244141
  cx a1,r
  ry r,-0.01017381500244141
  cx c4,a1
  cx c4,r
  ry r,0.01017381500244141
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.01017381500244141
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 1.7578125, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.10616154785156254
  cx a1,r
  ry r,-0.10616154785156254
  cx c3,a1
  cx c3,r
  ry r,0.10616154785156254
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-0.10616154785156254
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 0.12078825000000004 and theta = 0.1171875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,0.007077436523437503
  cx a2,r
  ry r,-0.007077436523437503
  cx c4,a2
  cx c4,r
  ry r,0.007077436523437503
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-0.007077436523437503
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.ry_c_x_to_the_p
  ##### Performs Ry(c*(x^p)) for c=0.08699629910714292 and p=7
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 1.0, controlled by [0]
  ry r,0.04349814955357146
  cx c0,r
  ry r,-0.08699629910714292
  cx c0,r
  ry r,0.04349814955357146
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.0078125, controlled by [1]
  ry r,0.000339829293387277
  cx c1,r
  ry r,-0.000679658586774554
  cx c1,r
  ry r,0.000339829293387277
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 6.103515625e-05, controlled by [2]
  ry r,2.6549163545881018e-06
  cx c2,r
  ry r,-5.3098327091762035e-06
  cx c2,r
  ry r,2.6549163545881018e-06
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 4.76837158203125e-07, controlled by [3]
  ry r,2.0741534020219545e-08
  cx c3,r
  ry r,-4.148306804043909e-08
  cx c3,r
  ry r,2.0741534020219545e-08
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 3.725290298461914e-09, controlled by [4]
  ry r,1.620432345329652e-10
  cx c4,r
  ry r,-3.240864690659304e-10
  cx c4,r
  ry r,1.620432345329652e-10
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 16.078125, controlled by [0, 1]
  ry r,0.6993686857910161
  cx c1,r
  ry r,-0.6993686857910161
  cx c0,c1
  cx c0,r
  ry r,0.6993686857910161
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-0.6993686857910161
  cx c0,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.1256103515625, controlled by [1, 2]
  ry r,0.005463817857742313
  cx c2,r
  ry r,-0.005463817857742313
  cx c1,c2
  cx c1,r
  ry r,0.005463817857742313
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.005463817857742313
  cx c1,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.029439926147460938, controlled by [1, 3]
  ry r,0.0012805823104083548
  cx c3,r
  ry r,-0.0012805823104083548
  cx c1,c3
  cx c1,r
  ry r,0.0012805823104083548
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.0012805823104083548
  cx c1,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.010005444288253784, controlled by [1, 4]
  ry r,0.00043521831200039044
  cx c4,r
  ry r,-0.00043521831200039044
  cx c1,c4
  cx c1,r
  ry r,0.00043521831200039044
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.00043521831200039044
  cx c1,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.00022999942302703857, controlled by [2, 4]
  ry r,1.0004549300065272e-05
  cx c4,r
  ry r,-1.0004549300065272e-05
  cx c2,c4
  cx c2,r
  ry r,1.0004549300065272e-05
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-1.0004549300065272e-05
  cx c2,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.0009813308715820312, controlled by [2, 3]
  ry r,4.268607701361182e-05
  cx c3,r
  ry r,-4.268607701361182e-05
  cx c2,c3
  cx c2,r
  ry r,4.268607701361182e-05
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-4.268607701361182e-05
  cx c2,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.5286306738853455, controlled by [0, 4]
  ry r,0.02299445611127002
  cx c4,r
  ry r,-0.02299445611127002
  cx c0,c4
  cx c0,r
  ry r,0.02299445611127002
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.02299445611127002
  cx c0,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 1.2806968688964844, controlled by [0, 3]
  ry r,0.055707943936049976
  cx c3,r
  ry r,-0.055707943936049976
  cx c0,c3
  cx c0,r
  ry r,0.055707943936049976
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.055707943936049976
  cx c0,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 7.666647434234619e-06, controlled by [3, 4]
  ry r,3.3348497666884235e-07
  cx c4,r
  ry r,-3.3348497666884235e-07
  cx c3,c4
  cx c3,r
  ry r,3.3348497666884235e-07
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-3.3348497666884235e-07
  cx c3,r
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 3.768310546875, controlled by [0, 2]
  ry r,0.1639145357322694
  cx c2,r
  ry r,-0.1639145357322694
  cx c0,c2
  cx c0,r
  ry r,0.1639145357322694
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.1639145357322694
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.5205960273742676, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.02264496385572107
  cx a0,r
  ry r,-0.02264496385572107
  cx c4,a0
  cx c4,r
  ry r,0.02264496385572107
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.02264496385572107
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 29.28515625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,1.2738501062622078
  cx a0,r
  ry r,-1.2738501062622078
  cx c2,a0
  cx c2,r
  ry r,1.2738501062622078
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-1.2738501062622078
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 3.242156982421875, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.14102782929754268
  cx a0,r
  ry r,-0.14102782929754268
  cx c3,a0
  cx c3,r
  ry r,0.14102782929754268
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.14102782929754268
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 11.524749755859375, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,0.5013052884478573
  cx a0,r
  ry r,-0.5013052884478573
  cx c3,a0
  cx c3,r
  ry r,0.5013052884478573
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.5013052884478573
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.228790283203125, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.009951953955173499
  cx a0,r
  ry r,-0.009951953955173499
  cx c3,a0
  cx c3,r
  ry r,0.009951953955173499
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.009951953955173499
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.0253293514251709, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.0011017799163870522
  cx a0,r
  ry r,-0.0011017799163870522
  cx c4,a0
  cx c4,r
  ry r,0.0011017799163870522
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0011017799163870522
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 5.112793922424316, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.22239707467420414
  cx a0,r
  ry r,-0.22239707467420414
  cx c4,a0
  cx c4,r
  ry r,0.22239707467420414
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.22239707467420414
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.09003710746765137, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.003916447565998885
  cx a0,r
  ry r,-0.003916447565998885
  cx c4,a0
  cx c4,r
  ry r,0.003916447565998885
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.003916447565998885
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.001787424087524414, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,7.774964027479296e-05
  cx a0,r
  ry r,-7.774964027479296e-05
  cx c4,a0
  cx c4,r
  ry r,7.774964027479296e-05
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-7.774964027479296e-05
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 1.4123454093933105, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.06143441183909033
  cx a0,r
  ry r,-0.06143441183909033
  cx c4,a0
  cx c4,r
  ry r,0.06143441183909033
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.06143441183909033
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.11640787124633789, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.005063526992686096
  cx a1,r
  ry r,-0.005063526992686096
  cx c4,a1
  cx c4,r
  ry r,0.005063526992686096
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.005063526992686096
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 2.8496646881103516, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,0.12395514078095564
  cx a1,r
  ry r,-0.12395514078095564
  cx c4,a1
  cx c4,r
  ry r,0.12395514078095564
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.12395514078095564
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 6.841678619384766, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.2976003597834708
  cx a1,r
  ry r,-0.2976003597834708
  cx c4,a1
  cx c4,r
  ry r,0.2976003597834708
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.2976003597834708
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 0.9281587600708008, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,0.04037318855501714
  cx a1,r
  ry r,-0.04037318855501714
  cx c4,a1
  cx c4,r
  ry r,0.04037318855501714
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.04037318855501714
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 14.90020751953125, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,0.6481314550638203
  cx a1,r
  ry r,-0.6481314550638203
  cx c3,a1
  cx c3,r
  ry r,0.6481314550638203
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-0.6481314550638203
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 0.08699629910714292 and theta = 2.58270263671875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,0.11234278554439553
  cx a2,r
  ry r,-0.11234278554439553
  cx c4,a2
  cx c4,r
  ry r,0.11234278554439553
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-0.11234278554439553
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.ry_c_x_to_the_p
  ##### Performs Ry(c*(x^p)) for c=0.07163903575086811 and p=9
  
  ### Rotations controlled by 1 qubit(s)
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 1.0, controlled by [0]
  ry r,0.035819517875434054
  cx c0,r
  ry r,-0.07163903575086811
  cx c0,r
  ry r,0.035819517875434054
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.001953125, controlled by [1]
  ry r,6.995999585045714e-05
  cx c1,r
  ry r,-0.00013991999170091427
  cx c1,r
  ry r,6.995999585045714e-05
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 3.814697265625e-06, controlled by [2]
  ry r,1.366406168954241e-07
  cx c2,r
  ry r,-2.732812337908482e-07
  cx c2,r
  ry r,1.366406168954241e-07
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 7.450580596923828e-09, controlled by [3]
  ry r,2.668762048738752e-10
  cx c3,r
  ry r,-5.337524097477504e-10
  cx c3,r
  ry r,2.668762048738752e-10
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 1.4551915228366852e-11, controlled by [4]
  ry r,5.212425876442875e-13
  cx c4,r
  ry r,-1.042485175288575e-12
  cx c4,r
  ry r,5.212425876442875e-13
  
  ### Rotations controlled by 2 qubit(s)
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 37.44140625, controlled by [0, 1]
  ry r,1.3411331204532633
  cx c1,r
  ry r,-1.3411331204532633
  cx c0,c1
  cx c0,r
  ry r,1.3411331204532633
  cx c1,r
  cx c0,c1
  cx c0,r
  ry r,-1.3411331204532633
  cx c0,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.07312774658203125, controlled by [1, 2]
  ry r,0.00261940062588528
  cx c2,r
  ry r,-0.00261940062588528
  cx c1,c2
  cx c1,r
  ry r,0.00261940062588528
  cx c2,r
  cx c1,c2
  cx c1,r
  ry r,-0.00261940062588528
  cx c1,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.012598782777786255, controlled by [1, 3]
  ry r,0.00045128232491762547
  cx c3,r
  ry r,-0.00045128232491762547
  cx c1,c3
  cx c1,r
  ry r,0.00045128232491762547
  cx c3,r
  cx c1,c3
  cx c1,r
  ry r,-0.00045128232491762547
  cx c1,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.003684585099108517, controlled by [1, 4]
  ry r,0.0001319800618210755
  cx c4,r
  ry r,-0.0001319800618210755
  cx c1,c4
  cx c1,r
  ry r,0.0001319800618210755
  cx c4,r
  cx c1,c4
  cx c1,r
  ry r,-0.0001319800618210755
  cx c1,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 2.460699761286378e-05, controlled by [2, 4]
  ry r,8.814107908547373e-07
  cx c4,r
  ry r,-8.814107908547373e-07
  cx c2,c4
  cx c2,r
  ry r,8.814107908547373e-07
  cx c4,r
  cx c2,c4
  cx c2,r
  ry r,-8.814107908547373e-07
  cx c2,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.00014282763004302979, controlled by [2, 3]
  ry r,5.116016847432187e-06
  cx c3,r
  ry r,-5.116016847432187e-06
  cx c2,c3
  cx c2,r
  ry r,5.116016847432187e-06
  cx c3,r
  cx c2,c3
  cx c2,r
  ry r,-5.116016847432187e-06
  cx c2,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.7256807258818299, controlled by [0, 4]
  ry r,0.025993533732582164
  cx c4,r
  ry r,-0.025993533732582164
  cx c0,c4
  cx c0,r
  ry r,0.025993533732582164
  cx c4,r
  cx c0,c4
  cx c0,r
  ry r,-0.025993533732582164
  cx c0,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 1.8865075707435608, controlled by [0, 3]
  ry r,0.06757379165239065
  cx c3,r
  ry r,-0.06757379165239065
  cx c0,c3
  cx c0,r
  ry r,0.06757379165239065
  cx c3,r
  cx c0,c3
  cx c0,r
  ry r,-0.06757379165239065
  cx c0,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 2.7896021492779255e-07, controlled by [3, 4]
  ry r,9.99222040514099e-09
  cx c4,r
  ry r,-9.99222040514099e-09
  cx c3,c4
  cx c3,r
  ry r,9.99222040514099e-09
  cx c4,r
  cx c3,c4
  cx c3,r
  ry r,-9.99222040514099e-09
  cx c3,r
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 6.4505767822265625, controlled by [0, 2]
  ry r,0.23105655035782424
  cx c2,r
  ry r,-0.23105655035782424
  cx c0,c2
  cx c0,r
  ry r,0.23105655035782424
  cx c2,r
  cx c0,c2
  cx c0,r
  ry r,-0.23105655035782424
  cx c0,r
  
  ### Rotations controlled by 3 qubit(s)
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 1.0835354402661324, controlled by [0, 3, 4]
  toffoli c0,c3,a0
  ry r,0.038811717071279034
  cx a0,r
  ry r,-0.038811717071279034
  cx c4,a0
  cx c4,r
  ry r,0.038811717071279034
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.038811717071279034
  cx c4,r
  toffoli c0,c3,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 108.9697265625, controlled by [0, 1, 2]
  toffoli c0,c1,a0
  ry r,3.90324306848663
  cx a0,r
  ry r,-3.90324306848663
  cx c2,a0
  cx c2,r
  ry r,3.90324306848663
  cx a0,r
  cx c2,a0
  cx c2,r
  ry r,-3.90324306848663
  cx c2,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 8.23084831237793, controlled by [0, 2, 3]
  toffoli c0,c2,a0
  ry r,0.29482501825520746
  cx a0,r
  ry r,-0.29482501825520746
  cx c3,a0
  cx c3,r
  ry r,0.29482501825520746
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.29482501825520746
  cx c3,r
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 38.66721153259277, controlled by [0, 1, 3]
  toffoli c0,c1,a0
  ry r,1.3850408746848966
  cx a0,r
  ry r,-1.3850408746848966
  cx c3,a0
  cx c3,r
  ry r,1.3850408746848966
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-1.3850408746848966
  cx c3,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.2128314971923828, controlled by [1, 2, 3]
  toffoli c1,c2,a0
  ry r,0.007623521618137949
  cx a0,r
  ry r,-0.007623521618137949
  cx c3,a0
  cx c3,r
  ry r,0.007623521618137949
  cx a0,r
  cx c3,a0
  cx c3,r
  ry r,-0.007623521618137949
  cx c3,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.016075875610113144, controlled by [1, 3, 4]
  toffoli c1,c3,a0
  ry r,0.0005758301137797021
  cx a0,r
  ry r,-0.0005758301137797021
  cx c4,a0
  cx c4,r
  ry r,0.0005758301137797021
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0005758301137797021
  cx c4,r
  toffoli c1,c3,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 16.338426545262337, controlled by [0, 1, 4]
  toffoli c0,c1,a0
  ry r,0.5852345616944905
  cx a0,r
  ry r,-0.5852345616944905
  cx c4,a0
  cx c4,r
  ry r,0.5852345616944905
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.5852345616944905
  cx c4,r
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.07552189752459526, controlled by [1, 2, 4]
  toffoli c1,c2,a0
  ry r,0.0027051579583689386
  cx a0,r
  ry r,-0.0027051579583689386
  cx c4,a0
  cx c4,r
  ry r,0.0027051579583689386
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.0027051579583689386
  cx c4,r
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.0004156865179538727, controlled by [2, 3, 4]
  toffoli c2,c3,a0
  ry r,1.4889690660425682e-05
  cx a0,r
  ry r,-1.4889690660425682e-05
  cx c4,a0
  cx c4,r
  ry r,1.4889690660425682e-05
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-1.4889690660425682e-05
  cx c4,r
  toffoli c2,c3,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 3.3820099756121635, controlled by [0, 2, 4]
  toffoli c0,c2,a0
  ry r,0.12114196677633618
  cx a0,r
  ry r,-0.12114196677633618
  cx c4,a0
  cx c4,r
  ry r,0.12114196677633618
  cx a0,r
  cx c4,a0
  cx c4,r
  ry r,-0.12114196677633618
  cx c4,r
  toffoli c0,c2,a0
  
  ### Rotations controlled by 4 qubit(s)
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 0.163043774664402, controlled by [1, 2, 3, 4]
  toffoli c1,c2,a0
  toffoli c3,a0,a1
  ry r,0.00584014940106979
  cx a1,r
  ry r,-0.00584014940106979
  cx c4,a1
  cx c4,r
  ry r,0.00584014940106979
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.00584014940106979
  cx c4,r
  toffoli c3,a0,a1
  toffoli c1,c2,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 13.789967447519302, controlled by [0, 1, 3, 4]
  toffoli c0,c1,a0
  toffoli c3,a0,a1
  ry r,0.49394998548807134
  cx a1,r
  ry r,-0.49394998548807134
  cx c4,a1
  cx c4,r
  ry r,0.49394998548807134
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.49394998548807134
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 36.644615828990936, controlled by [0, 1, 2, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,1.3125924717249544
  cx a1,r
  ry r,-1.3125924717249544
  cx c4,a1
  cx c4,r
  ry r,1.3125924717249544
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-1.3125924717249544
  cx c4,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 3.450474813580513, controlled by [0, 2, 3, 4]
  toffoli c0,c2,a0
  toffoli c3,a0,a1
  ry r,0.12359434426378217
  cx a1,r
  ry r,-0.12359434426378217
  cx c4,a1
  cx c4,r
  ry r,0.12359434426378217
  cx a1,r
  cx c4,a1
  cx c4,r
  ry r,-0.12359434426378217
  cx c4,r
  toffoli c3,a0,a1
  toffoli c0,c2,a0
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 83.47841262817383, controlled by [0, 1, 2, 3]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  ry r,2.9901564933477323
  cx a1,r
  ry r,-2.9901564933477323
  cx c3,a1
  cx c3,r
  ry r,2.9901564933477323
  cx a1,r
  cx c3,a1
  cx c3,r
  ry r,-2.9901564933477323
  cx c3,r
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  
  ### Rotations controlled by 5 qubit(s)
  
  # Rotation c*theta with c = 0.07163903575086811 and theta = 22.6483154296875, controlled by [0, 1, 2, 3, 4]
  toffoli c0,c1,a0
  toffoli c2,a0,a1
  toffoli c3,a1,a2
  ry r,0.8112517393821603
  cx a2,r
  ry r,-0.8112517393821603
  cx c4,a2
  cx c4,r
  ry r,0.8112517393821603
  cx a2,r
  cx c4,a2
  cx c4,r
  ry r,-0.8112517393821603
  cx c4,r
  toffoli c3,a1,a2
  toffoli c2,a0,a1
  toffoli c0,c1,a0
  

.result
  display