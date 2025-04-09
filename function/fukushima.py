import numpy as np
from tripledif import *
import math             # na binomický koeficient math.comb()

# Funkcia Fukushima (WIP), cieľom je výpočet gravitačného potenciálu
# zrýchlenia, a tenzoru na základe článku Fukushima 2018. 

# TO DO : 
#       vyriešiť indexy a opraviť výpočty koeficientov c vo výpočte V,g. zatiaľ sú sumácie dané manuálne 
# 
#       dopočítať tenzory
#
#       prerobiť funkcie z formátu lambda na formát def fx return..
#
# POZNÁMKY : 
#       indexy j používať bežne bez riešenia čo sú, pri konštantnom polynóme je N = 0 ! 

def Fukushima(coor_prism,coor_point,rho):

    # INPUTs: 
    #           "coor_prism"        - 3D karteziánske súradnice okrajov pravouhlej prizmy vo tvare "x1,x2,y1,y2,z1,z2" (vektor s 6 hodnotami)
    #           "coor_point"        - 3D karteziánske súradnice výpočtového bodu (vektor s 3 hodnotami)
    #           "rho"               - vektor parametrov polynómu vyjadrujúceho hustotu (zmenu) v radiálnom smere prizmy (rôzne dlhý vektor)

    # OUTPUTs:
    #                               Zatiaľ nič - do budúcnosti gravitačný potenciál, zrýchlenie a tenzor

    G = 6.67430 * 10**-11           # Newtonova gravitačná konštanta (zdr. Bucha, Fyzikálna Geodézia )     

    # ================================== F U N K C I A ================================================== #

    X1, Y1, Z1 = coor_prism[0] - coor_point[0], coor_prism[2] - coor_point[1], coor_prism[4] - coor_point[2]                # (10) posunutie "shifted endpoints"
    X2, Y2, Z2 = coor_prism[1] - coor_point[0], coor_prism[3] - coor_point[1], coor_prism[5] - coor_point[2] 
    
    N = len(rho) - 1       # stupeň polynómu

    # ---------- Homogénna prizma ------------- # 

    # Koeficienty polynómu -------

    c_mj = math.comb(0,0) * G * rho[0]                                                                                     # (15)
    d_c0_nm = (0 + 1) * math.comb((0 + 0 + 1) , 0) * G * rho[0]                                                            # (19)
    dd_c0_nm = (0 + 2) * (0 + 1) * math.comb((0 + 0 + 2), 0) * G * rho[0]

    c_0 = c_mj * (coor_point[2]**0)                                                                                        # (14)
    d_c_0 = d_c0_nm * (coor_point[2]**0)
    dd_c_0 = dd_c0_nm * (coor_point[2]**0)

    # Potenciál - V ------
    U_0 = lambda X, Y, Z, func_val_0: \
       -((X**2 * func_val_0[0]) + (Y**2 * func_val_0[1]) + (Z**2 * func_val_0[2])) / 2 + \
             (Y * Z * func_val_0[3]) + (Z * X * func_val_0[4]) + (X * Y * func_val_0[5])                                   # (23) ako anonýmna funkcia lambda

    # def U_0(X,Y,Z):
    #     return -((X**2 * elemfun_0[0] + Y**2 * elemfun_0[1] + Z**2 * elemfun_0[2])/2) + \
    #             Y * Z * elemfun_0[3] + Z * X * elemfun_0[4] + X * Y * elemfun_0[5]



    W_0 = triple_dif(U_0,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)                                                                 # (9)

    V = (c_0 * W_0)                                                                                                         # (13)

    # Zrýchlenie - g ------ 

    U_0x = lambda X, Y, Z, func_val_0: \
            -(X * func_val_0[0]) + (Y * func_val_0[5]) + (Z * func_val_0[4])                                                # (23)

    U_0y = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[5]) - (Y * func_val_0[1]) + (Z * func_val_0[3])
    
    U_0z = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[4]) + (Y * func_val_0[3]) - (Z * func_val_0[2])
    
    W_0X = -(triple_dif(U_0x,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))                                                           # (20)
    W_0Y = -(triple_dif(U_0y,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))
    W_0Z = -(triple_dif(U_0z,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))

    # Počiatočná hodnota zrýchlenia ----

    g_x = c_0 * W_0X
    g_y = c_0 * W_0Y
    g_z = (c_0 * W_0Z) # + (d_c_0 * W_0)  # zakomentované, Bucha to vo výpočte odstránil - opýtať sa prečo (asi kvôli sum^-1)

    g = np.array([g_x,g_y,g_z])   # vloženie do vektora (tuple)

    # Tenzor - T ------

    U_0XX = lambda X, Y, Z, func_val_0: \
                -func_val_0[0]
    
    U_0XY = lambda X, Y, Z, func_val_0: \
                func_val_0[5]
    
    U_0XZ = lambda X, Y, Z, func_val_0: \
                func_val_0[4]
    
    U_0YY = lambda X, Y, Z, func_val_0: \
                -func_val_0[1]
    
    U_0YZ = lambda X, Y, Z, func_val_0: \
                func_val_0[3]
    
    U_0ZZ = lambda X, Y, Z, func_val_0: \
                -func_val_0[2]
    
    W_0XX = triple_dif(U_0XX,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)
    W_0XY = triple_dif(U_0XY,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)
    W_0XZ = triple_dif(U_0XZ,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)
    W_0YY = triple_dif(U_0YY,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)
    W_0YZ = triple_dif(U_0YZ,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)
    W_0ZZ = triple_dif(U_0ZZ,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)

    # Počiatočná hodnota tenzora ----

    G_0XX = c_0 * W_0XX
    G_0XY = c_0 * W_0XY
    G_0XZ = (c_0 * W_0XZ) # + (d_c_0 * W_0X)
    G_0YY = c_0 * W_0YY
    G_0YZ = (c_0 * W_0YZ) # + (d_c_0 * W_0Y)
    G_0ZZ = (c_0 * W_0ZZ) # + (2 * d_c_0 * W_0Z) + (dd_c_0 * W_0) # ten istý dôvod zakomentovania

    G     = np.array([[G_0XX, G_0XY, G_0XZ],    # vloženie do matice
                      [0    , G_0YY, G_0YZ],
                      [0    , 0    , G_0ZZ]])

    # ---------- Nehomogénna prizma ------------- # 

    # :) 

     #   U_1 = lambda X, Y, Z, func_val_1: \
     #   -((Z**(1+2) * func_val_1[5]) / (1 + 2)) + ((Z**(1+1) * ((Y * func_val_1[2]) + (X * func_val_1[1]))) / (1 + 1)) \
     #   - (((Y * func_val_1[4]) + (X * func_val_1[3])) / ((1 + 1) * (1 + 2)))
     # 
     #   W_1 = triple_dif(U_1,X1,Y1,Z1,X2,Y2,Z2,N,homogenous=False)        

    return V, g, G




