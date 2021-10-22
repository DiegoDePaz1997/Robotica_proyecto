from defs_robotica import *
import numpy as np


#T_s(b1, a2, a3, d6, ang)
b1 = 1
a2 = 1
a3 = 1
d6 = .1
ang= (1.1, 0.2, -1.3, 0.4, 1.5, 0.6,)

T = T_s(b1, a2, a3, d6, ang)
tt = T[5]
print(tt)


tt= np.array([
[0.2487,0.8039,-0.5403,-0.0540],
[0.1597,0.5162,0.8415,0.0841],
[0.9553,-0.2955,0.0000,0.0000],
[0.0000,0.0000,0.0000,1.0000]
])
angulo = cin_inversa(b1, a2, a3, d6, tt)
ang= (1.1, .2, -1.3   , .4, 1.5, .6)
T = T_s(b1, a2, a3, d6, ang)


ang456 = []
for i in range(4):
    t = angulo[i]
   
    A1 =  DH(0, np.pi/2,  b1, t[0])
    A2 =  DH(a2,     0 ,  0, t[1])
    A3 =  DH(a3,     0,   0, t[2])
    T03 = A1@A2@A3
    T03 = T03.transpose()
    np.set_printoptions(suppress=True, precision=4)
    ang456.append(angulos(T03@tt))
   


ang456 = np.round(ang456, 4)



final =[list(angulo[0]) + list(ang456[0][0] ),
list(angulo[0]) + list(ang456[0][1] ),
list(angulo[1]) + list(ang456[1][0] ),
list(angulo[1]) + list(ang456[1][1] ),

list(angulo[2]) + list(ang456[2][0] ),
list(angulo[2]) + list(ang456[2][1] ),
list(angulo[3]) + list(ang456[3][1] ),
list(angulo[3]) + list(ang456[3][0] )
]

[print(i) for i in final]





T = T_s(b1, a2, a3, d6, final[0])
tt = T[5]
print(tt)


T = T_s(b1, a2, a3, d6, final[4])
tt = T[5]
print(tt)



