import numpy as np

# Doplňujúca knižnica tripledif ku knižnici fukushima

# zoznam článkov (zdroje vzťahov) je uvedený v knižnici fukushima

# Aktuálne obsahuje 3 funkcie - 
#                               elemfun                - elementárne funkcie (n = 0,1,2..N)
#                               triple_dif             - operátor trojitej diferencie (Fukushima (11))

# TO DO: 
#            skontrolovať správnosť (obzvlášť pri elemfun_n  - rekurentný výpočet)
#            upraviť triple_dif aby obsahoval argument stupňa polynómu N bez potreby *args

def atan3(X,Y,Z):
    # Funkcia na výpočet špeciálnej funkcie atan3 podľa (37) kapitoly 2.6 (Fukushima)

    if X == 0:
        return 0
    else:
       return np.arctan((Y * Z) / (X * np.sqrt(X**2 + Y**2 + Z**2 )))  
    
def logsum(X,Y,Z):
    # Funkcia na výpočet špeciálnej funkcie logsum podľa (38) kapitoly 2.6 (Fukushima)

    omega = 1e-150                                          # (39) -    konštanta omega
    if X > 0:
        return np.log(X + np.sqrt(X**2 + Y**2 + Z**2))
    elif X == 0:
        return (1/2) * np.log(Y**2 + Z**2 + omega)
    else:
        return np.log((Y**2 + Z**2 + omega) / (np.sqrt(X**2 + Y**2 + Z**2) - X))

def elemfun(X,Y,Z):

    # funkcia pre výpočet elementárnych funkcií A,B,C,D,E,F (24) v prípade homogénnej prizmy 

    # INPUT -   "X,Y,Z"                 - trojrozmerné súradnice bodu posunuté podľa (10) 

    # OUTPUT -  "A,B,C,D,E,F"           - elementárne funkcie podľa (24) v prípade homogénnej prizmy  
    #           "R"                     - euklidovská vzdialenosť

    # ============================ F U N K C I A =================================== #

    omega = 1e-150                                          # (39) -    konštanta omega

    A = atan3(X,Y,Z)
    B = atan3(Y,Z,X)
    C = atan3(Z,X,Y)
    
    D = logsum(X,Y,Z)
    E = logsum(Y,Z,X)
    F = logsum(Z,X,Y)
    R = np.sqrt(X**2 + Y**2 + Z**2 + omega)                  # (25) -    euklidovská vzdialenosť

    return A, B, C, D, E, F, R
    #      0, 1, 2, 3, 4, 5, 6

def triple_dif(equation,X1,Y1,Z1,X2,Y2,Z2):

    # funkcia trojitej diferencie aplikuje operátor podľa vzťahu Fukushima (11) 

    # INPUT -   "equation"              - rovnica obsahujúca súradnice X,Y,Z zadefinovaná anonýmnou funkciou lambda, viď. DEFINE
    #           "X1,Y1...Z2"            - upravené súradnice, shifted endpoints (10)
    #           "*args"                 - dodatočné argumenty
    #           "N"                     - stupeň polynómu modelujúceho hustotu, hodnota 0 uvažuje s homogénnou prizmou

    # OUTPUT -  "triple_difference"     - vyčíslená hodnota zvolenej rovnice s použitím operátora trojitej diferencie 

    # DEFINE -  "equation"              - rovnica musí byť definovaná cez lambda (anonýmnu funkciu), kde parametrami funkcie sú súradnice XYZ a 
    #                                      tuple hodnoty vypočítaných element. funkcií pre danú prizmu v poradí v akom ich funkcia elemfun vracia 

    # ============================ F U N K C I A =================================== #

    f1_element_val = elemfun(X2,Y2,Z2)          # výpočet elementárnych funkcií pre rovnicu
    f1 = equation(X2, Y2, Z2,f1_element_val)    # výpočet rovnice na základe daných súradníc 

    f2_element_val = elemfun(X2,Y2,Z1)
    f2 = equation(X2, Y2, Z1,f2_element_val)

    f3_element_val = elemfun(X2,Y1,Z2)
    f3 = equation(X2, Y1, Z2,f3_element_val)

    f4_element_val = elemfun(X2,Y1,Z1)
    f4 = equation(X2, Y1, Z1,f4_element_val)

    f5_element_val = elemfun(X1,Y2,Z2)
    f5 = equation(X1, Y2, Z2,f5_element_val)

    f6_element_val = elemfun(X1,Y2,Z1)
    f6 = equation(X1, Y2, Z1,f6_element_val)

    f7_element_val = elemfun(X1,Y1,Z2)
    f7 = equation(X1, Y1, Z2,f7_element_val)

    f8_element_val = elemfun(X1,Y1,Z1)
    f8 = equation(X1, Y1, Z1,f8_element_val)


    triple_difference = f1 - f2 - f3 + f4 - f5 + f6 + f7 - f8            # (11) -    trojitá diferencia funkcie (operátor)

    return triple_difference