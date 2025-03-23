import numpy as np
from tripledif import *
import math             # na binomický koeficient math.comb()

# Funkcia Fukushima (WIP), cieľom je výpočet gravitačného potenciálu
# zrýchlenia, a tenzoru na základe článku Fukushima 2018. 

# TO DO : 
#
# vyriešiť indexy a opraviť výpočty koeficientov c vo výpočte V,g. zatiaľ sú sumácie dané manuálne 
# 
# dopočítať tenzory

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
    
    N = len(rho)        # stupeň polynómu

    # ---------- Homogénna prizma ------------- # 

    # Potenciál - V ------
    U_0 = lambda X, Y, Z, func_val_0: \
        -((X**2 * func_val_0[0]) + (Y**2 * func_val_0[1]) + (Z**2 * func_val_0[2]) / 2) + \
              (Y * Z * func_val_0[3]) + (Z * X * func_val_0[4]) + (X * Y * func_val_0[5])                                   # (23) ako anonýmna funkcia lambda
    U_1 = lambda X, Y, Z, func_val_1: \
        -((Z**(1+2) * func_val_1[5]) / (1 + 2)) + ((Z**(1+1) * ((Y * func_val_1[2]) + (X * func_val_1[1]))) / (1 + 1)) \
        - (((Y * func_val_1[4]) + (X * func_val_1[3])) / ((1 + 1) * (1 + 2)))                                               # (26)

    W_0 = triple_dif(U_0,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)                                                                 # (9)
    W_1 = triple_dif(U_1,X1,Y1,Z1,X2,Y2,Z2,N,homogenous=False)

    # V tomto bode je výpočet samotných hodnôt potenciálu a zrýchlenia veľmi pochybný pretože si nie som istý indexami,
    # z tohto titulu je výpočet neistý, všetky parametre a spôsob výpočtu parametrov c treba prerobiť (aj pre zrýchlenie g)
    # až budú indexy isté 

    c_00 = math.comb(0,0) * G * rho[0]                                                                                      # (15)
    c_01 = math.comb(1,0) * G * rho[0]
    c_10 = math.comb(0,1) * G * rho[0]
    c_11 = math.comb(1,1) * G * rho[0]

    c_0 = (c_00 * (coor_point[2]**0)) + (c_01 * (coor_point[2]**1))                                                         # (14)
    c_1 = (c_10 * (coor_point[2]**0)) + (c_11 * (coor_point[2]**1))

    V = (c_0 * W_0) + (c_1 * W_1)                                                                                           # (13)

    # Zrýchlenie - g ------ 

    U_0x = lambda X, Y, Z, func_val_0: \
            -(X * func_val_0[0]) + (Y * func_val_0[5]) + (Z * func_val_0[4])                                                # (23)

    U_0y = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[5]) - (Y * func_val_0[1]) + (Z * func_val_0[3])
    
    U_0z = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[4]) + (Y * func_val_0[3]) - (Z * func_val_0[2])
    
    dW_0x = -(triple_dif(U_0x,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))                                                           # (20)
    dW_0y = -(triple_dif(U_0y,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))
    dW_0z = -(triple_dif(U_0z,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))

    # Tu je to pravdepodobne chybné pretože stále nerozumiem indexom, prerobiť až sa v nich vyznáš 
    # zatiaľ uvažujem že n = 1, m = 0 - ale úprimne neviem úplne prečo, just a hunch, aj pri potenciáli to mám inak :) 
    d_c0_nm = 1 * math.comb(1 + 0 + 1,1) * G * rho[0]                                                                       # (19)

    # Opäť som odignoroval sumáciu kvôli indexom, zároveň indexy samotné :) 
    d_c_0 = d_c0_nm * coor_point[2]**1

    # Počiatočná hodnota zrýchlenia ----

    g_x = c_0 * dW_0x
    g_y = c_0 * dW_0y
    g_z = (c_0 * dW_0z) + (d_c_0 * W_0)

    # Tenzor - T ------


    

    # ---------- Nehomogénna prizma ------------- # 

    # :) 

    # -- Posledné úpravy --

    g = [g_x,g_y,g_z]   # vloženie do vektora (tuple)

    return V, g




