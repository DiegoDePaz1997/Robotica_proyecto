import numpy as np


#inversa de una matriz 4x4 como traspuesta
def inv_4x4(m):
    m = np.array(m)
    m[:3, :3] = m[:3, :3].transpose()
    m[:3, 3:]= -m[:3, :3]@m[:3, 3:]

    return m


if(__name__=='__main__'):   
    mat = [[1, 2,  3,  4], 
          [ 1, 2,  3,  5],
          [ 1, 2,  3,  6],
          [ 0, 0,  0,  1]]

    b1 = 1
    a2 = 1
    a3 = 1
    d6 = .1
    s = inv_4x4(mat)
 
    print(s)


#el @ es para multiplicacion matricial




#radianes_o_grados = 'radianes'
def radianes_o_grados(x):
    global grados
    grados = x


def dec_rad_o_deg_Rn(f):
    def inner(a):
        if(grados == 'radianes'): pass 
        else: a = np.radians(a)
        return f(a)
    return inner




'''
matrices puras de rotacion y traslacion
'''

def traslacion(d):
    # d = dx, dy, dz
    matriz = np.array([ [1, 0, 0, d[0]],
                        [0, 1, 0, d[1]],
                        [0, 0, 1, d[2]],
                        [0, 0, 0,   1]])
    return matriz    



@dec_rad_o_deg_Rn
def Rx(ang):
    matriz = np.array([[1, 0, 0, 0],
                        [0, np.cos(ang), -np.sin(ang), 0],
                        [0, np.sin(ang), np.cos(ang), 0],
                        [0, 0, 0, 1]])
    return matriz


@dec_rad_o_deg_Rn
def Ry(ang):
    matriz = np.array([[np.cos(ang), 0, np.sin(ang), 0],
                        [0, 1, 0, 0],
                        [-np.sin(ang), 0, np.cos(ang), 0],
                        [0, 0, 0, 1]])
    return matriz


@dec_rad_o_deg_Rn
def Rz(ang):
    matriz = np.array([[np.cos(ang), -np.sin(ang), 0, 0],
                        [np.sin(ang), np.cos(ang), 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])
    return matriz



'''
angulos de euler
'''
'''
Dada una matriz devuelve las 2 soluciones de los ángulos
theta phi y psi posibles
'''
def angulos(matriz):  #ZYZ
    #devuelve dos soluciones en forma de matriz
    cos_theta = matriz[2][2]    
    sen_theta = (1 - cos_theta**2)**(1/2)    
    cos_psi, sin_psi = matriz[0, 0], matriz[0, 1]    

    if(abs(sen_theta)< 10**-6):

        if(cos_theta > 0):
            theta_1 = 0
            phi_1 = 0
            psi_1 = np.arctan2(-sin_psi, cos_psi)
        else:
            theta_1 = np.pi
            phi_1 = 0
            psi_1 = -np.arctan2(-sin_psi, -cos_psi)

        theta_2 = theta_1  
        phi_2 = phi_1
        psi_2 = psi_1

        
    else:  
        phi_1   = np.arctan2(matriz[1, 2], matriz[0, 2])
        theta_1 =  np.arctan2(sen_theta, cos_theta)
        psi_1   = np.arctan2(matriz[2, 1], -matriz[2,0])

        phi_2   = np.arctan2(-matriz[1, 2], -matriz[0, 2])
        theta_2 =  np.arctan2(-sen_theta, cos_theta)
        psi_2   = np.arctan2(-matriz[2, 1], matriz[2, 0])

    a = [phi_1, theta_1, psi_1]
    b = [phi_2, theta_2, psi_2]
    return [a, b]

lista = [
[0.2487,	0.8039,	-0.5403,	-0.0540],
[0.1597,	0.5162,	0.8415,	0.0841],
[0.9553,	-0.2955,	0,	0],
[0,	0,	0,	1]
]
lista = np.array(lista)
print(angulos(lista))

#####################################################
def R_ZYZ(angulos):
    return Rz(angulos[0])@(Ry(angulos[1])@Rz(angulos[2]))    
#####################################################


##########################################################################################################



'''
Representacion de eje/angulo
'''

#dado un vector k y un ángulo encontrar la matriz de rotación resultante
def R_K_phi(k, phi):
    #k tiene que ser un vector unitario
    kx, ky, kz = k
    vo = 1 - np.cos(phi)
    matriz = np.array([
        [kx**2*vo + np.cos(phi)   , kx*ky*vo - kz*np.sin(phi), kx*kz*vo + ky*np.sin(phi)],
        [kx*ky*vo + kz*np.sin(phi), ky**2*vo + np.cos(phi)   , ky*kz*vo - kx*np.sin(phi)],
        [kx*kz*vo - ky*np.sin(phi), ky*kz*vo + kx*np.sin(phi), kz**2*vo +np.cos(phi)]
    ])
    return matriz
    


'''
dada una matriz encontrar phi y K
'''
def K_phi_dada_R(arr):
    phi = np.arccos((np.trace(arr) - 1)/2)# -2 era -1 la ponerla de 4x4 la puse en -2
    
    matriz = np.array([arr[2][1] - arr[1][2],
                       arr[0][2] - arr[2][0],
                       arr[1][0] - arr[0][1]])
    
    K = matriz/(2*np.sin(phi))

    #dos soluciones k phi y -k -phi
    return K, phi
    
'''
dados alfa y betha encontrar la matriz
'''
def R_Alf_Bet(alpha, betha, phi):
    #phi rotacion del sistema
    R_k = Rz(alpha)@Ry(betha)@Rz(phi)@Ry(-betha)@Rz(-alpha)
    return R_k







'''
funcion denavit hartenberg
Cinematica directa
'''

def Tx(alpha, a): # a = eje x
    matriz = np.array([ [1, 0, 0, a],
                        [0, np.cos(alpha), -np.sin(alpha), 0],
                        [0, np.sin(alpha), np.cos(alpha), 0],
                        [0, 0, 0, 1]])
    return matriz    


def Tz(theta, b): # b = eje z
    matriz = np.array([[np.cos(theta), -np.sin(theta), 0, 0],
                        [np.sin(theta), np.cos(theta), 0, 0],
                        [0, 0, 1, b],
                        [0, 0, 0, 1]])
    return matriz


def DH(a, alpha, b, theta):
    radianes_o_grados('radianes')
    #en dos matrices juta lo que en la primera se logra en cuatro
    return Tz(theta, b)@Tx(alpha, a)
    #return traslacion((0,0,b))@Rz(theta)@Rx(alpha)@traslacion((a,0,0))



'''
DH aplicado a un robot
'''
def T_s(b1, a2, a3, d6, ang):

    # lista de parametros de DH
    tupla = (
        (0, np.pi/2, b1 , ang[0]),    # A1
        (a2, 0      , 0 , ang[1]),  # A2
        (a3, 0      , 0 , ang[2]),   # A3

        (0, -np.pi/2, 0 , ang[3]),     # A4
        (0,  np.pi/2, 0 , ang[4]),     # A5
        (0,  0      , d6, ang[5])     # A6
    )
    #genera las matrices 
    Ai = [DH(a, alpha, b, theta) for a, alpha, b, theta in tupla]


    aux = np.identity(4)
    lista_aux = []
    for Ti in Ai:
        aux = aux@Ti
        lista_aux.append(aux)
        '''
        T1 = A1
        T2 = A1*A2
        T3 = A1*A2*A3
        T4 = A1*A2*A3*A4
        T5 = A1*A2*A3*A4*A5
        T6 = A1*A2*A3*A4*A5*A6
        '''        

    return np.array(lista_aux)



'''
cinematica inversa 
la primera parte se resuelve con angulos de euler
'''
def cin_inversa(b1, a2, a3, d6, matriz):
    #angulos(matriz)
    xyz1 = matriz@np.array([0, 0, -d6, 1])
    Xc, Yc, Zc = xyz1[:-1] # posicion final del brazo 3

    if (Xc**2 + Yc**2)**.5 > a2+a3:
         print('error de alcance')
         return 0
    if abs(Xc)<10**-4 and abs(Yc)<10**-4:
         print('Xc y Yc no estan definidos en el origen')
         return 0


    theta1_pos = np.arctan2(Yc, Xc)

    theta1_neg = np.arctan2(-Yc, -Xc)

    c3 = (Xc**2 +Yc**2 + (Zc - b1)**2 - a2**2 - a3**2)/(2*a2*a3)
    s3 = (1-c3**2)**.5 #dos soluciones

    theta3_pos = np.arctan2(s3,  c3)
    theta3_neg = np.arctan2(-s3, c3)



    beta = np.arctan2(Zc-b1, (Xc**2 + Yc**2)**.5)

    gama_pos = np.arctan2( a3*np.sin(theta3_pos),
                           a2 + a3*c3)

    gama_neg =  np.arctan2( a3*np.sin(theta3_neg),
                            a2 + a3*c3)

    theta2_pos = beta - gama_pos
    theta2_neg = beta - gama_neg

    """otros valor de theta1"""


    lis = (
        #sol theta1
        (theta1_pos,theta2_pos, theta3_pos),
         (theta1_pos,theta2_neg, theta3_neg),

         #sol theta_pi
         (theta1_neg,theta2_pos, theta3_pos),
         (theta1_neg,theta2_neg, theta3_neg)
    )


    ang456 = []
    for i in range(4):
        t = lis[i]
    
        A1 =  DH(0, np.pi/2,  b1, t[0])
        A2 =  DH(a2,     0 ,  0, t[1])
        A3 =  DH(a3,     0,   0, t[2])
        T03 = inv_4x4(A1@A2@A3)
        ang456.append(angulos(T03@matriz))
    


    array = np.array(
    [ 
            [theta1_pos, theta2_pos, theta3_pos] + ang456[0][0],
            [theta1_pos, theta2_neg, theta3_neg] + ang456[1][0],
            [theta1_neg, theta2_pos, theta3_pos] + ang456[2][0],
            [theta1_neg, theta2_neg, theta3_neg] + ang456[3][0],

            [theta1_pos, theta2_pos, theta3_pos] + ang456[0][1],
            [theta1_pos, theta2_neg, theta3_neg] + ang456[1][1],
            [theta1_neg, theta2_pos, theta3_pos] + ang456[2][1],
            [theta1_neg, theta2_neg, theta3_neg] + ang456[3][1]

    ]      
    )
    array = array.round(4)




    return array










  
   
