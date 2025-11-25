import numpy as np

# Doplňujúca knižnica tripledif ku knižnici fukushima

# zoznam článkov (zdroje vzťahov) je uvedený v knižnici fukushima

# Aktuálne obsahuje 3 funkcie - 
#                               elemfun_0              - elementárne funkcie (n = 0)
#                               elemfun_n              - elementárne funkcie (n = 1,2..)
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
        return np.log((Y**2 + Z**2 + omega) / np.sqrt(X**2 + Y**2 + Z**2 - X))

def elemfun_0(X,Y,Z,N):

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

    # D = np.log(X + R)                                       # (24) -    elementárne funkcie
    # E = np.log(Y + R)
    # F = np.log(Z + R)
    # 
    # A = np.arctan((Y*Z) / (X*R))
    # B = np.arctan((Z*X) / (Y*R))
    # C = np.arctan((X*Y) / (Z*R))

    return A, B, C, D, E, F, R
    #      0, 1, 2, 3, 4, 5, 6

def elemfun_n(X,Y,Z,N):   

    # funkcia pre výpočet elementárnych funkcií v prípade nehomogénnej prizmy teda n = 1,2.. 

    # INPUT -   "X,Y,Z"                 - trojrozmerné súradnice bodu posunuté podľa (10)
    #           "N"                     - maximálna hodnota rekurentného výpočtu (stupeň polynómu)

    # OUTPUT -  "R_n, E_n, D_n"         - elementárne funkcie vypočítané pre daný stupeň polynómu
    #           "E_n2, D_n2"            - elementárne funkcie vypočítané pre stupeň polynómu + 2 pre potrebu výpočtu U 
    #           "C"                     - elementárna funkcia podľa (24) ktorá je taktiež potrebná pre výpočet U ak n=1,2.. 

    # ============================ F U N K C I A =================================== #
    A, B, C, D, E, F, R = elemfun_0(X, Y, Z,N)            # uloženie premenných z funkcie elemfun_0 pre výpočet v tejto funkcii 

    S = X**2 + Y**2                                                                 # (31)

    R_1 = R                                                                         # (30)
    R_2 = ((Z * R) - (S * F))/2
    
    D_1 = D                                                                         # (28)
    D_2 = (Y * B) - (X * F)

    E_1 = E
    E_2 = (X * A) - (Y * F)

    # Rekurentný výpočet R,D,E =====================================================
    if N == 1:       # v prípade že N = 1 ---------------------------------------
        
        D_n2 = -(Y**2 * D_1) - (X * R_1)                                             # (27)
        E_n2 = -(X**2 * E_1) - (Y * R_1)

        R_n = R_1       # ukladanie premenných pre potrebu výpočtu U
        D_n = D_1
        E_n = E_1


    elif N == 2:     # v prípade že N = 2 ---------------------------------------

        D_n2 = -(Y**2 * D_2) - (X * R_2)                                             # (27)
        E_n2 = -(X**2 * E_2) - (Y * R_2)

        R_n = R_2
        D_n = D_2
        E_n = E_2

    else:               # v prípade že N > 2 ---------------------------------------
        for n in range(3,N + 3):     # range od 3 (0,1,2 máme) po nmax + 3 (+1 pretože python neberie poslednú hodnotu, + 2 pretože potrebuje n+2 na výpočet U)
            
            R_n = ((Z**(n-1) * R) - ((n - 1) * S * R_1)) / n                        # (29) 
            D_n = -(Y**2 * D_1) - (X * R_1)                                         # (27)
            E_n = -(X**2 * E_1) - (Y * R_1)

            if n == (N+1):
                R_main = R_n
                D_main = D_n
                E_main = E_n

            R_1 = R_2       # prepísanie premenných 
            R_2 = R_n 

            D_1 = D_2
            D_2 = D_n 

            E_1 = E_2
            E_2 = E_n 
        # koniec cyklu for - - - - - 
    
        R_n = R_main        # ukladanie premenných pre potrebu výpočtu U

        D_n2 = D_n 
        D_n = D_main

        E_n2 = E_n 
        E_n = E_main

    return R_n, E_n, D_n, E_n2, D_n2, C 
    #      0    1    2    3     4     5



def triple_dif(equation,X1,Y1,Z1,X2,Y2,Z2,N,*args):

    # funkcia trojitej diferencie aplikuje operátor podľa vzťahu Fukushima (11) 

    # INPUT -   "equation"              - rovnica obsahujúca súradnice X,Y,Z zadefinovaná anonýmnou funkciou lambda, viď. DEFINE
    #           "X1,Y1...Z2"            - upravené súradnice, shifted endpoints (10)
    #           "*args"                 - dodatočné argumenty
    #           "N"                     - stupeň polynómu modelujúceho hustotu, hodnota 0 uvažuje s homogénnou prizmou

    # OUTPUT -  "triple_difference"     - vyčíslená hodnota zvolenej rovnice s použitím operátora trojitej diferencie 

    # DEFINE -  "equation"              - rovnica musí byť definovaná cez lambda (anonýmnu funkciu), kde parametrami funkcie sú súradnice XYZ a 
    #                                      tuple hodnoty vypočítaných element. funkcií pre danú prizmu v poradí v akom ich funkcia elemfun vracia 

    # ============================ F U N K C I A =================================== #

    if N == 0:                      # Predefinovanie názvu funkcie, použitie element. funkc. pre homogénnu či nehomogénnu prizmu 
        elemfun = elemfun_0
    else:
        elemfun = elemfun_n

    f1_element_val = elemfun(X2,Y2,Z2,N,*args)    # výpočet elementárnych funkcií pre rovnicu
    f1 = equation(X2, Y2, Z2,f1_element_val)    # výpočet rovnice na základe daných súradníc 

    f2_element_val = elemfun(X2,Y2,Z1,N,*args)
    f2 = equation(X2, Y2, Z1,f2_element_val)

    f3_element_val = elemfun(X2,Y1,Z2,N,*args)
    f3 = equation(X2, Y1, Z2,f3_element_val)

    f4_element_val = elemfun(X2,Y1,Z1,N,*args)
    f4 = equation(X2, Y1, Z1,f4_element_val)

    f5_element_val = elemfun(X1,Y2,Z2,N,*args)
    f5 = equation(X1, Y2, Z2,f5_element_val)

    f6_element_val = elemfun(X1,Y2,Z1,N,*args)
    f6 = equation(X1, Y2, Z1,f6_element_val)

    f7_element_val = elemfun(X1,Y1,Z2,N,*args)
    f7 = equation(X1, Y1, Z2,f7_element_val)

    f8_element_val = elemfun(X1,Y1,Z1,N,*args)
    f8 = equation(X1, Y1, Z1,f8_element_val)


    triple_difference = f1 - f2 - f3 + f4 - f5 + f6 + f7 - f8            # (11) -    trojitá diferencia funkcie (operátor)

    return triple_difference