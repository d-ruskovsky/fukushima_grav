import numpy as np
from tripledif import *
import math             # na binomický koeficient math.comb()

# Funkcia Fukushima (WIP), cieľom je výpočet gravitačného potenciálu
# zrýchlenia, a tenzoru na základe článku Fukushima 2018. 

# TO DO : 
#
# Začať riešiť POTENCIAL, výpočet zatiaľ riešiť pre homogénnu prizmu (n=0), 
# výpočet aj tak musí začať homogénnou prizmou, až potom riešiť výpočet cez 
# cyklus pre rho = 2+. 
#
# Vo vzorci (13) vstupuje súradnica z výpočtového bodu, vo vzorci (14) je táto
# súradnica umocnená na index j ktorý = n - m
#
# OTÁZNE - index j netreba riešiť, stačí ho nahradiť n - m. N je stupeň 
# polynómu rho, ??? m nadobúda hodnoty 0-N a tiež n ??? 
#
# Bolo by dobré túto funkciu prepísať do jazyka c, ale to až bude funkčná

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
    
    # ---------- Homogénna prizma ------------- # 

    # Potenciál - V ------
    U_0 = lambda X, Y, Z, func_val: \
        -((X**2 * func_val[0]) + (Y**2 * func_val[1]) + (Z**2 * func_val[2]) / 2) + \
              (Y * Z * func_val[3]) + (Z * X * func_val[4]) + (X * Y * func_val[5])                                         # (23) ako anonýmna funkcia lambda

    W_0 = triple_dif(U_0,X1,Y1,Z1,X2,Y2,Z2,homogenous=True)                                                                 # (9)

    # Indexy !! prekontrolovať tieto vzorce až budú indexy vyriešené
    c_0j = 1 * G * rho[0]                                                                                                   # (15) pri n=0,  binom. koef = 1 

    c_0 = c_0j * (coor_point[2]**1)                                                                                         # (14)

    # Počiatočná hodnota potenciálu ----
    V = c_0 * W_0                                                                                                           # (13)

    # Zrýchlenie - g ------ 

    U_0x = lambda X, Y, Z, func_val: \
            -(X * func_val[0]) + (Y * func_val[5]) + (Z * func_val[4])                                                      # (23)

    U_0y = lambda X, Y, Z, func_val: \
            (X * func_val[5]) - (Y * func_val[1]) + (Z * func_val[3])
    
    U_0z = lambda X, Y, Z, func_val: \
            (X * func_val[4]) + (Y * func_val[3]) - (Z * func_val[2])
    
    dW_0x = -(triple_dif(U_0x,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))                                                           # (20)
    dW_0y = -(triple_dif(U_0y,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))
    dW_0z = -(triple_dif(U_0z,X1,Y1,Z1,X2,Y2,Z2,homogenous=True))

    # Tu je to pravdepodobne chybné pretože stále nerozumiem indexom, prerobiť až sa v nich vyznáš 
    # zatiaľ uvažujem že n = 1, m = 0 - ale úprimne neviem úplne prečo, just a hunch, aj pri potenciáli to mám inak :) 
    d_c0_nm = 1 * math.comb(1 + 0 + 1,1) * G * rho[0]                                                                       # (19)

    # Opäť som odignoroval sumáciu kvôli indexom, zároveň indexy :) 
    d_c_0 = d_c0_nm * coor_point[2]**1

    # Počiatočná hodnota zrýchlenia ----

    g_x = c_0 * dW_0x
    g_y = c_0 * dW_0y
    g_z = (c_0 * dW_0z) + (d_c_0 * W_0)

    # ---------- Nehomogénna prizma ------------- # 

    # :) 

    # -- Posledné úpravy --

    g = [g_x,g_y,g_z]

    return V, g




