import numpy as np
'''
angulos de euler
'''


'''
defs Rx Ry Rz
'''
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


@dec_rad_o_deg_Rn
def Rx(ang):
    matriz = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])
    return matriz


@dec_rad_o_deg_Rn
def Ry(ang):
    matriz = np.array([[np.cos(ang), 0, np.sin(ang)],
                        [0, 1, 0],
                        [-np.sin(ang), 0, np.cos(ang)]])
    return matriz


@dec_rad_o_deg_Rn
def Rz(ang):
    matriz = np.array([[np.cos(ang), -np.sin(ang), 0],
                        [np.sin(ang), np.cos(ang), 0],
                        [0, 0, 1]])
    return matriz




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



#####################################################
def R_ZYZ(angulos):
    return Rz(angulos[0])@(Ry(angulos[1])@Rz(angulos[2]))    
#####################################################


##########################################################################################################



'''
Representacion de eje/angulo
'''

#dado un vector k y un ángulo encontrar la matriz de rotación resultante
def R_K_phi(kx, ky, kz, phi):
    vo = 1 - np.cos(phi)
    matriz = np.array([
        [kx**2*vo + np.cos(phi)   , kx*ky*vo - kz*np.sin(phi), kx*kz*vo + ky*np.sin(phi)],
        [kx*ky*vo + kz*np.sin(phi), ky**2*vo + np.cos(phi)   , ky*kz*vo - kx*np.sin(phi)],
        [kx*kz*vo - ky*np.sin(phi), ky*kz*vo + kx*np.sin(phi), kz**2*vo +np.cos(phi)]
    ])

    


'''
dada una matriz encontrar phi y K
'''
def K_phi_dada_R(arr):
    phi = np.arccos((np.trace(arr) - 1)/2)
    
    matriz = np.array([[arr[2][1] - arr[1][2]],
                  [arr[0][2] - arr[2][0]],
                  [arr[1][0] - arr[0][1]]])
    
    K = matriz/(2*np.sin(phi))

    #dos soluciones k phi y -k -phi
    return (K, phi)
    













'''
    convierte radianes a grados
    alfa = np.degrees(np.pi)
    beta = np.degrees(np.pi)
    theta = np.degrees(np.pi)

    convierte grados a radianes
    alfa = np.radians(60)
    beta = np.radians(30)
    theta = np.radians(90)
'''


#el @ es para multiplicacion matricial

if(__name__=='__main__'):   
    punto = np.array([1, 2, 3]) 

    radianes_o_grados('grados') 
    phi = 90
    theta = 90
    psi = 90   
    R_ = R_ZYZ((phi, theta, psi))
    
    
    print('angulos reales')
    print([phi, theta, psi])
    print(' ')
    ###################################
    print('angulos predichos')
    A = angulos(R_)
    print(A[0])# solucion 1
    print(A[1])# solucion 2
    print(' ')
    #diferencia entre las 2 soluciones casi cero
    print(R_ZYZ(A[0]) - R_ZYZ(A[1])) 


    radianes_o_grados('grados')
    print(' ')
    r = Ry(20)@Rz(30)@Rx(25)@Ry(15)
    print(r)

    z = K_phi_dada_R(r)
    print(z)











  
   
