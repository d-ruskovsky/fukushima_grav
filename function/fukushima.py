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
#       skúsiť zjednotiť indexy m,n  pri c_coef - alebo lepšie pochopiť indexy
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

    # Elementárne funkcie U sú transformované na W pomocou funkcie triple_dif, ktorá si volá funkcie elemfun
    # Keďže tieto funkcie vyžadujú rôzne kombinácie parametrov a súradníc je na ich fungovanie potrebné 
    # zadefinovať elementárne funkcie U pomocou anonýmnej funkcie lambda alebo pomocou bežnej funkcie def, 
    # musia však mať definované parametre pomocou func_val v poradí v akom elemfun vracia funkcie ABCDEFR.
    # Nižšie sú preto tieto parametre zoradené, pre lepšie pochopenie zostavenia vzťahov U v samotnom výpočte. 

    #  elemfun
    #  A, B, C, D, E, F, R
    #  0, 1, 2, 3, 4, 5, 6

    # ================================== F U N K C I A ================================================== #

    X1, Y1, Z1 = coor_prism[0] - coor_point[0], coor_prism[2] - coor_point[1], coor_prism[4] - coor_point[2]                # (10) posunutie "shifted endpoints"
    X2, Y2, Z2 = coor_prism[1] - coor_point[0], coor_prism[3] - coor_point[1], coor_prism[5] - coor_point[2] 

    # ---------- Homogénna prizma ------------- # 

    # Váhové funkcie  - počiatočné hodnoty

    # Potenciál - V ------
    # U = lambda X, Y, Z, func_val_0: \
    #    -((X**2 * func_val_0[0]) + (Y**2 * func_val_0[1]) + (Z**2 * func_val_0[2])) / 2 + \
    #          (Y * Z * func_val_0[3]) + (Z * X * func_val_0[4]) + (X * Y * func_val_0[5])                     # (23) ako anonýmna funkcia lambda
    
    def U(X, Y, Z, func_val):                                                                                  # (23) funckia prerobená 
        return(
            -((X**2 * func_val[0]) + (Y**2 * func_val[1]) + (Z**2 * func_val[2])) / 2 +
            (Y * Z * func_val[3]) + (Z * X * func_val[4]) + (X * Y * func_val[5])
        )

    # Zrýchlenie - g ------ 
    Ux = lambda X, Y, Z, func_val_0: \
            -(X * func_val_0[0]) + (Y * func_val_0[5]) + (Z * func_val_0[4])                                 # (23)

    Uy = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[5]) - (Y * func_val_0[1]) + (Z * func_val_0[3])
    
    Uz = lambda X, Y, Z, func_val_0: \
            (X * func_val_0[4]) + (Y * func_val_0[3]) - (Z * func_val_0[2])
    
    # Tenzor - T ------
    Uxx = lambda X, Y, Z, func_val_0: \
                -func_val_0[0]
    
    Uxy = lambda X, Y, Z, func_val_0: \
                func_val_0[5]
    
    Uxz = lambda X, Y, Z, func_val_0: \
                func_val_0[4]
    
    Uyy = lambda X, Y, Z, func_val_0: \
                -func_val_0[1]
    
    Uyz = lambda X, Y, Z, func_val_0: \
                func_val_0[3]
    
    Uzz = lambda X, Y, Z, func_val_0: \
                -func_val_0[2]
    
    # Inicializácia gravitačného účinku do cyklu
    V,g_x,g_y,g_z,G_XX,G_XY,G_YY,G_XZ,G_YZ,G_ZZ  = (0,) * 10

    for q in range(N+1):
        # Kontrola homogenity - ak je q 0 tak je prizma homogénna - použijú sa vzťahy pre homogénnu prizmu
        if q > 0:
            
            def U(X, Y, Z, func_val):
                return (
                    -((Z**(q+2) * func_val[2]) / (q+2))
                    + ((Z**(q+1) * (Y*func_val[3] + X*func_val[4])) / (q+1))
                    - ((Y*func_val[3] + X*func_val[4]) / ((q+1)*(q+2))) )
            
            Ux = lambda X, Y, Z, func_val: \
            (Z**(q+1)*func_val[4] - func_val[4]) / (q+1)

            Uy = lambda X, Y, Z, func_val: \
            (Z**(q+1)*func_val[3] - func_val[3]) / (q+1)

            Uz = lambda X, Y, Z, func_val: \
            -Z**(q+1) * func_val[2] + Z**q * (Y * func_val[3] + X * func_val[4])

            Uxx = lambda X, Y, Z, func_val: \
            X * func_val[4]

            Uxy = lambda X, Y, Z, func_val: \
            func_val[6]

            Uxz = lambda X, Y, Z, func_val: \
            Z**q * func_val[4]

            Uyy = lambda X, Y, Z, func_val: \
            Y * func_val[3]

            Uyz = lambda X, Y, Z, func_val: \
            Z**q * func_val[3]

            Uzz = lambda X, Y, Z, func_val: \
            -(q + 1)*Z**q * func_val[2] + q * Z**(q-1) * (Y*func_val[3] + X*func_val[4])
        
        # Aplikácia operátora trojitej diferencie
        
        W = triple_dif(U,X1,Y1,Z1,X2,Y2,Z2)                                                                # (9)                                                                                                    
       
        Wx = -(triple_dif(Ux,X1,Y1,Z1,X2,Y2,Z2))                                                           # (20)
        Wy = -(triple_dif(Uy,X1,Y1,Z1,X2,Y2,Z2))
        Wz = -(triple_dif(Uz,X1,Y1,Z1,X2,Y2,Z2))  

        Wxx = triple_dif(Uxx,X1,Y1,Z1,X2,Y2,Z2)
        Wxy = triple_dif(Uxy,X1,Y1,Z1,X2,Y2,Z2)
        Wxz = triple_dif(Uxz,X1,Y1,Z1,X2,Y2,Z2)
        Wyy = triple_dif(Uyy,X1,Y1,Z1,X2,Y2,Z2)
        Wyz = triple_dif(Uyz,X1,Y1,Z1,X2,Y2,Z2)
        Wzz = triple_dif(Uzz,X1,Y1,Z1,X2,Y2,Z2)

        # Výpočet koeficientov polynómu
        c, dc, ddc = c_coef(N,q,q,Z2,rho_s,rho_gradient)

        # Výpočet hodnôt gravitačného účinku (pripočítanie ku iniciálnej hodnote)
        V += c * W

        g_x += c * Wx
        g_y += c * Wy
        g_z += (c * Wz) + (dc * W)

        G_XX += c * Wxx
        G_XY += c * Wxy
        G_YY += c * Wyy
        G_XZ += (c * Wxz) + (dc * Wx)
        G_YZ += (c * Wyz) + (dc * Wy)
        G_ZZ += (c * Wzz) + (dc * Wz) + (ddc * W)

    # Vloženie do array
    g = np.array([g_x,g_y,g_z])
    
    G = np.array([[G_XX, G_XY, G_XZ],
                  [0    , G_YY, G_YZ], 
                  [0    , 0    , G_ZZ]])
        
    return V, g, G

def c_coef(N, m, n, z, rho_s, rho_gradient):
    # INPUTs: 
    #           "N"                 - stupeň polynómu modelujúceho hustotu
    #           "m"                 - aktuálny krok pri sumácií výpočtu grav. potenciálu
    #           "n"                 - aktuálny krok pri sumácií výpočtu g rav. zrýchlenia a tenzora
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
    for j in range(0, N - m + 1):
        rho = density(rho_s, rho_gradient, z, j+m)
        cmj = math.comb(j + m, m) * G * rho                                     # (15)        # !!! používam z ako depth, skonzultovať či je to správne - takmer určite nie
        c += cmj * (z**j)                                                       # (14)

    # Koeficient c'  - v tomto bode sú v článku koeficienty označované ako c'nm c''nm, ja ich označujem ako c'nj c''nj keďže m používa koef. c, vzťahy sú inak nemenné - taktiež asi nesprávne 
    if N - n - 1 >= 0:
        for j in range(0,N - n - 1 + 1):
            rho = density(rho_s, rho_gradient, z, n + j + 1)
            dcnj = (j + 1) * math.comb(n + j + 1, n) * G * rho                                           # (19)
            dc += dcnj * (z**j)                                                                          # (18)
    else:
        dc = 0


    # Koeficient c''
    if N - n - 2 >= 0:
        for j in range(0,N - n - 2 + 1):
            rho = density(rho_s, rho_gradient, z, n + j + 2)
            ddcnj = (j + 2) * (j + 1) * math.comb(n + j + 2, n) * G * rho                                # (19)
            ddc += ddcnj * (z**j)                                                                       # (18)
    else:
        ddc = 0


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