import numpy as np
from tripledif import * # dotatočná knižnica 
import math             # na binomický koeficient math.comb()

# Knižnica fukushima (WIP) 
# cieľom je výpočet gravitačného účinku (gravitačný potenciál, zrýchlenie, tenzor) pravouhlej prizmy s vertikálne premenlivou hustotou
 
# Vzťahy obsiahnuté v tomto skripte sú prevzaté z článkov:
#
#   1.   Fukushima 2018 : Recursive computation of gravitational field of a right rectangular parallelepiped with density varying vertically by following an arbitrary degree polynomial
#   2.   Goossens et al 2020 : High-resolution gravity field models from GRAIL data and implications for models of the density structure of the Moons crust

# Knižnica aktuálne obsahuje funkcie:
#       fukushima()     - hlavná funkcia, výpočet gravitačného účinku prizmy na výpočtový bod
#       c_coef()        - výpočet koeficientov polynómu hustoty
#       density()       - výpočet hustoty pre daný polynóm

# TO DO LIST: 
# 
#       upraviť tripledif - nepáči sa mi potreba písania N= , a ani args*
#
#       prerobiť funkcie z formátu lambda na formát def fx return..
#
#       skúsiť zjednotiť indexy m,n  pri c_coef
#       
#       prekonzultovať rozvinutie hustoty nad stupeň 1
#
#       prekonzultovať rozdiel medzi súradnicou z   a hĺbkou depth

def Fukushima(coor_prism,coor_point,rho_s,rho_gradient,N):

    # INPUTs: 
    #           "coor_prism"        - 3D karteziánske súradnice okrajov pravouhlej prizmy vo tvare "x1,x2,y1,y2,z1,z2" (vektor s 6 hodnotami), resp. doména prizmy
    #           "coor_point"        - 3D karteziánske súradnice výpočtového bodu (vektor s 3 hodnotami)
    #           "rho_s"             - hustota na povrchu prizmy v jednotách kg/m3
    #           "rho_gradient"      - gradient hustoty v jednotkách kg/m3/m
    #           "N"                 - stupeň polynómu modelujúceho hustotu

    # OUTPUTs:
    #           "V"                 - gravitačný potenciál generovaný prizmou vo výpočtovom bode
    #           "g"                 - vektor gravitačného zrýchlenia generujúceho prizmou vo výpočtovom bode
    #           "G"                 - tenzor gravitačného zrýchlenia generujúceho prizmou vo výpočtovom bode  

    # ================================== F U N K C I A ================================================== #

    X1, Y1, Z1 = coor_prism[0] - coor_point[0], coor_prism[2] - coor_point[1], coor_prism[4] - coor_point[2]                # (10) posunutie "shifted endpoints"
    X2, Y2, Z2 = coor_prism[1] - coor_point[0], coor_prism[3] - coor_point[1], coor_prism[5] - coor_point[2] 

    # ---------- Homogénna prizma ------------- # 

    # Váhové funkcie

    # Potenciál - V ------
    U_0 = lambda X, Y, Z, func_val_0: \
       -((X**2 * func_val_0[0]) + (Y**2 * func_val_0[1]) + (Z**2 * func_val_0[2])) / 2 + \
             (Y * Z * func_val_0[3]) + (Z * X * func_val_0[4]) + (X * Y * func_val_0[5])                     # (23) ako anonýmna funkcia lambda

    # Zrýchlenie - g ------ 
    U_0x = lambda X, Y, Z, func_val_0: \
            -(X * func_val_0[0]) + (Y * func_val_0[5]) + (Z * func_val_0[4])                                 # (23)

    U_0y = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[5]) - (Y * func_val_0[1]) + (Z * func_val_0[3])
    
    U_0z = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[4]) + (Y * func_val_0[3]) - (Z * func_val_0[2])
    
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
    
    # Aplikácia operátora trojitej diferencie

    W_0 = triple_dif(U_0,X1,Y1,Z1,X2,Y2,Z2,N=0)                                                                # (9)                                                                                                    # (13)

    W_0X = -(triple_dif(U_0x,X1,Y1,Z1,X2,Y2,Z2,N=0))                                                           # (20)
    W_0Y = -(triple_dif(U_0y,X1,Y1,Z1,X2,Y2,Z2,N=0))
    W_0Z = -(triple_dif(U_0z,X1,Y1,Z1,X2,Y2,Z2,N=0))  
    
    W_0XX = triple_dif(U_0XX,X1,Y1,Z1,X2,Y2,Z2,N=0)
    W_0XY = triple_dif(U_0XY,X1,Y1,Z1,X2,Y2,Z2,N=0)
    W_0XZ = triple_dif(U_0XZ,X1,Y1,Z1,X2,Y2,Z2,N=0)
    W_0YY = triple_dif(U_0YY,X1,Y1,Z1,X2,Y2,Z2,N=0)
    W_0YZ = triple_dif(U_0YZ,X1,Y1,Z1,X2,Y2,Z2,N=0)
    W_0ZZ = triple_dif(U_0ZZ,X1,Y1,Z1,X2,Y2,Z2,N=0)

    # Výpočet koeficientov polynómu
    c0, dc0, ddc0 = c_coef(N,0,0,Z2,rho_s,rho_gradient) # v tomto prípade by mali byť koeficienty dc0 a ddc0 nulové

    # Výpočet počiatočných hodnôt gravitačného účinku (homogénna prizma)
    V0 = c0 * W_0

    g_x0 = c0 * W_0X
    g_y0 = c0 * W_0Y
    g_z0 = (c0 * W_0Z) + (dc0 * W_0)

    G_XX0 = c0 * W_0XX
    G_XY0 = c0 * W_0XY
    G_YY0 = c0 * W_0YY
    G_XZ0 = (c0 * W_0XZ) + (dc0 * W_0X)
    G_YZ0 = (c0 * W_0YZ) + (dc0 * W_0Y)
    G_ZZ0 = (c0 * W_0ZZ) + (dc0 * W_0Z) + (ddc0 * W_0)

    if N == 0:          # v prípade homogénnej prizmy vrátiť hodnoty pre N=0 a ukončiť program
        V = V0

        g = np.array([g_x0,g_y0,g_z0])

        G = np.array([[G_XX0, G_XY0, G_XZ0],
                      [0    , G_YY0, G_YZ0], 
                      [0    , 0    , G_ZZ0]])
        
        # print("Calculating homogenous prism")
        return V, g, G
    else:
        # ---------- Nehomogénna prizma ------------- # 
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Tu si skončil naposledy, treba doprogramovať skript na výpočet g. účinku nehomogénnej prizmy, 
        # ostatné funkcie by mali fungovať fajn, pred tým ako začneš, zamysli sa nad indexmi v cm, konkrétne 
        # N - m, hoc pri N=1 môžeš použiť W_0 pri sumácií, nemôžeš použiť c0 pretože N už nie je 0 ale 1. 
        # nájdi preto nejaký rozumný spôsob ako pripojiť W_0 k sumácií nehomogénnej prizmy (c_coef musia 
        # byť použité nanovo), zároveň to skús future-proofnuť ak by si v budúcnosti bol moc premotivovaný 
        # a chcel skúsiť exponenciálnu hustotu, alebo vyšší polynóm abo niečo také.


     #   U_1 = lambda X, Y, Z, func_val_1: \
     #   -((Z**(1+2) * func_val_1[5]) / (1 + 2)) + ((Z**(1+1) * ((Y * func_val_1[2]) + (X * func_val_1[1]))) / (1 + 1)) \
     #   - (((Y * func_val_1[4]) + (X * func_val_1[3])) / ((1 + 1) * (1 + 2)))
     # 
     #   W_1 = triple_dif(U_1,X1,Y1,Z1,X2,Y2,Z2,N,homogenous=False)     

     # g = np.array([g_x,g_y,g_z])   # vloženie do vektora (tuple)
     # 
     # G = np.array([[G_xx, G_xy, G_xz],     # xx   xy    xz       # vloženie do matice
     #               [0   , G_yy, G_yz],     # yx   yy    yz
     #               [0   , 0   , G_zz]])    # zx   zy    zz

        return None


def c_coef(N, m, n, z, rho_s, rho_gradient):
    # INPUTs: 
    #           "N"                 - stupeň polynómu modelujúceho hustotu
    #           "m"                 - aktuálny krok pri sumácií výpočtu grav. potenciálu
    #           "n"                 - aktuálny krok pri sumácií výpočtu grav. zrýchlenia a tenzora
    #           "z"                 - súradnice z prizmy (aktuálne taktiež hĺbka - depth)
    #           "rho_s"             - hustota na povrchu prizmy v jednotkách kg/m3
    #           "rho_gradient"      - gradient hustoty v jednotkách kg/m3/m

    # OUTPUTs:
    #           "c, dc, ddc"        - koeficient (a jeho derivácie) polynómu hustoty 
    
    G = 6.67430 * 10**-11           # Newtonova gravitačná konštanta (zdr. Bucha, Fyzikálna Geodézia 2023)   

    c = 0
    dc = 0
    ddc = 0

    # Koeficient c
    for j in range(N - m + 1):
        rho = density(rho_s, rho_gradient, z, j+m)
        cmj = math.comb(j + m, m) * G * rho                                     # (15)        # !!! používam z ako depth, skonzultovať či je to správne 
        c = c + (cmj * z**j)                                                    # (14)

    # Koeficient c'  - v tomto bode sú v článku koeficienty označované ako c'nm c''nm, ja ich označujem ako c'nj c''nj keďže m používa koef. c, vzťahy sú inak nemenné
    if (N - n - 1) < 0:
        dc = 0
    else:
        for j in range(N - n - 1 + 1):
            rho = density(rho_s, rho_gradient, z, n + j + 1)
            dcnj = (j + 1) * math.comb(n + j + 1, n) * G * rho                                           # (19)
            dc = dc + (dcnj * z**j)                                                                      # (18)

    # Koeficient c''
    if (N - n - 2) < 0:
        ddc = 0
    else:
        for j in range(N - n - 2 + 1):
            rho = density(rho_s, rho_gradient, z, n + j + 2)
            ddcnj = (j + 2) * (j + 1) * math.comb(n + j + 2, n) * G * rho                                # (19)
            ddc = ddc + (ddcnj * z**j)                                                                   # (18)

    return c, dc, ddc


def density(rho_s, rho_gradient, depth, N):
    # INPUTs: 
    #           "rho_s"             - hustota na povrchu prizmy v jednotách kg/m3
    #           "rho_gradient"      - gradient hustoty v jednotkách kg/m3/m
    #           "depth"             - hĺbka, resp. obrátená súradnice z počítanej prizmy
    #           "N"                 - stupeň polynómu modelujúceho hustotu

    # OUTPUTs:
    #           "rho"               - výsledná hustota pre daný stupeň, hĺbku, gradient a hustotu na povrchu
    
    # TO DO:
    #
    # - konzultovať a dokončiť polynóm pri požadovanom stupni vyššom ako 1 


    if N == 0:                                  # prípad homogénnej prizmy - vráti hustotu na povrchu
        return rho_s
    
    elif N == 1:                                # prípad prizmy s lineárne premenlivou hustotou
        rho = rho_s + rho_gradient * depth      # (6) - Goossens et al 2020
        return rho
    
    else:                                       # PLACEHOLDER - v prípade potreby rho vyššieho stupňa ako 1 sa vráti stupeň 1
        rho = rho_s + rho_gradient * depth      # tento krok treba prekonzultovať (rozšíriť na exponenciálnu funkciu? 
        return rho                              # alebo rozšíriť Taylorom? prípadne nechať tento stupeň ?