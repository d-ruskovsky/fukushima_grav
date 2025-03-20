import numpy as np
from tripledif import *

# Funkcia Fukushima (WIP), cieľom je výpočet gravitačného potenciálu
# zrýchlenia, a tenzoru na základe článku Fukushima 2018. 

# TO DO : 
#
# Začať riešiť POTENCIAL, výpočet zatiaľ riešiť pre homogénnu prizmu (n=0), 
# výpočet aj tak musí začať homogénnou prizmou, až potom riešiť výpočet cez 
# cyklus pre rho = 2+. 
#
# Vytvoriť funkciu pre riešenie triple difference (ten trojuholnik 3) ? 
#
# Vo vzorci (13) vstupuje súradnica z výpočtového bodu, vo vzorci (14) je táto
# súradnica umocnená na index j ktorý = n - m
#
# OTÁZNE - index j netreba riešiť, stačí ho nahradiť n - m. N je stupeň 
# polynómu rho, ??? m nadobúda hodnoty 0-N a tiež n ??? 
#
# Bolo by dobré túto funkciu prepísať do jazyka c, ale to až bude funkčná
#
# ??????!!!!!!!!! Netreba ako vstup použiť skôr krajné "súradnice" prizmy x1,x2,ako dx ??? do výpočtu vstupujú viac a nikde nie je žiadny rozmer prizmy

def Fukushima(x1,y1,z1,dx, dy, dz, poloha_bod, poloha_element,rho):

    # INPUTs: 
    #           "poloha_bod"        - 3D karteziánske súradnice výpočtového bodu
    #                               (vektor s 3 hodnotami)
    #           "poloha_element"    - poloha elementu (?) v doméne prizmy, v článku označené ako x',y',z'
    #           "rho"               - vektor parametrov polynómu vyjadrujúceho 
    #                               hustotu (zmenu) v radiálnom smere prizmy
    #           "dx,dy,dz"          - rozmery prizmy v 3D kart. súr. sys. 
    #           "x1,x2,x3"          - súradnice okrajov prizmy ? asi zle, opýtať sa !!!!

    # OUTPUTs:
    #                               Zatiaľ nič - do budúcnosti gravitačný potenciál, zrýchlenie a tenzor
    G = 6.67430 * 10**-11           # Newtonova gravitačná konštanta (zdr. Bucha, Fyzikálna Geodézia )     

    # ================================== F U N K C I A ==================================================
    
    x2, y2, z2 = x1 + dx, y1 + dy, z1 + dz                                                                                  # pripočítanie súradníc rohov prizmy podľa jej zadanej veľkosti (?)

    X, Y, Z = poloha_element[0] - poloha_bod[0], poloha_element[1] - poloha_bod[1], poloha_element[2] - poloha_bod[2]       # (6) zmena integračných premenných (asi netreba)

    X1, Y1, Z1 = x1 - poloha_bod[0], y1 - poloha_bod[1], z1 - poloha_bod[2]                                                 # (10) posunutie "shifted endpoints"
    X2, Y2, Z2 = x2 - poloha_bod[0], y2 - poloha_bod[1], z2 - poloha_bod[2] 

    U_0 = "-((X**2 * A) + (Y**2 * B) + (Z**2 * C) / 2) + (Y * Z * D) + (Z * X * E) + (X * Y * F)"                           # (23) v stringu pre funkciu triple_dif

    W_0 = triple_dif(U_0,X1,Y1,Z1,X2,Y2,Z2)                                                                                 # (9)

    c_0j = 1 * G * rho[0]                                                                                                   # (15) tieto neviem či sú správne stále nechápem indexom m,n,j
                                                                                                                            # ale tipujem že bin. koef. je pri konštantnom polynome rho 1 

    c_0 = c_0j * (poloha_bod[2]**1)                                                                                         # (14)


    c = c_0
    W = W_0
    V = c * W                                                                                                               # (13)

    return V




