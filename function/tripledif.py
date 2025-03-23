import numpy as np

# Knižnica tripledif 

# Aktuálne obsahuje 3 funkcie - 
#                               elemfun_0              - elementárne funkcie (n = 0)
#                               elemfun_n              - elementárne funkcie (n = 1,2..)
#                               triple_dif             - operátor trojitej diferencie (Fukushima (11))

# TO DO: 
#            skontrolovať správnosť (obzvlášť pri elemfun_n  - rekurentný výpočet)
#            upraviť triple_dif aby obsahoval argument stupňa polynómu N bez potreby *args

def elemfun_0(X,Y,Z):

    # funkcia pre výpočet elementárnych funkcií A,B,C,D,E,F (24) v prípade homogénnej prizmy 

    # INPUT -   "X,Y,Z"                 - trojrozmerné súradnice bodu posunuté podľa (10) 

    # OUTPUT -  "A,B,C,D,E,F"           - elementárne funkcie podľa (24) v prípade homogénnej prizmy  
    #           "R"                     - euklidovská vzdialenosť

    # ============================ F U N K C I A =================================== #

    R = np.sqrt(X**2 + Y**2 + Z**2)                         # (25) -    euklidovská vzdialenosť

    D = np.log(X + R)                                       # (24) -    elementárne funkcie
    E = np.log(Y + R)
    F = np.log(Z + R)
    
    A = np.arctan((Y*Z) / (X*R))
    B = np.arctan((Z*X) / (Y*R))
    C = np.arctan((X*Y) / (Z*R))

    return A, B, C, D, E, F, R
    #      0, 1, 2, 3, 4, 5, 6

def elemfun_n(X,Y,Z,nmax):   

    # funkcia pre výpočet elementárnych funkcií v prípade nehomogénnej prizmy teda n = 1,2.. 

    # INPUT -   "X,Y,Z"                 - trojrozmerné súradnice bodu posunuté podľa (10)
    #           "nmax"                  - maximálna hodnota rekurentného výpočtu (stupeň polynómu)

    # OUTPUT -  "R_n, E_n, D_n"         - elementárne funkcie vypočítané pre daný stupeň polynómu
    #           "E_n2, D_n2"            - elementárne funkcie vypočítané pre stupeň polynómu + 2 pre potrebu výpočtu U 
    #           "C"                     - elementárna funkcia podľa (24) ktorá je taktiež potrebná pre výpočet U ak n=1,2.. 

    # ============================ F U N K C I A =================================== #

    R = elemfun_0(X,Y,Z)[6]         # uloženie premenných z funkcie elemfun_0 pre výpočet v tejto funkcii 
    F = elemfun_0(X,Y,Z)[5]
    E = elemfun_0(X,Y,Z)[4]
    D = elemfun_0(X,Y,Z)[3]
    C = elemfun_0(X,Y,Z)[2]
    B = elemfun_0(X,Y,Z)[1]
    A = elemfun_0(X,Y,Z)[0]

    S = X**2 + Y**2                                                                 # (31)

    R_1 = R                                                                         # (30)
    R_2 = ((Z * R) - (S * F))/2
    
    D_1 = D                                                                         # (28)
    D_2 = (Y * B) - (X * F)

    E_1 = E
    E_2 = (X * A) - (Y * F)

    # Rekurentný výpočet R,D,E =====================================================
    if nmax == 1:       # v prípade že N = 1 ---------------------------------------
        
        D_n2 = -(Y**2 * D_1) - (X * R_1)                                             # (27)
        E_n2 = -(X**2 * E_1) - (Y * R_1)

        R_n = R_1       # ukladanie premenných pre potrebu výpočtu U
        D_n = D_1
        E_n = E_1


    elif nmax == 2:     # v prípade že N = 2 ---------------------------------------

        D_n2 = -(Y**2 * D_2) - (X * R_2)                                             # (27)
        E_n2 = -(X**2 * E_2) - (Y * R_2)

        R_n = R_2
        D_n = D_2
        E_n = E_2

    else:               # v prípade že N > 2 ---------------------------------------
        for n in range(3,nmax + 3):     # range od 3 (0,1,2 máme) po nmax + 3 (+1 pretože python neberie poslednú hodnotu, + 2 pretože potrebuje n+2 na výpočet U)
            
            R_n = ((Z**(n-1) * R) - ((n - 1) * S * R_1)) / n                        # (29) 
            D_n = -(Y**2 * D_1) - (X * R_1)                                         # (27)
            E_n = -(X**2 * E_1) - (Y * R_1)

            if n == (nmax+1):
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



def triple_dif(equation,X1,Y1,Z1,X2,Y2,Z2,*args,homogenous=True):

    # funkcia trojitej diferencie aplikuje operátor podľa vzťahu Fukushima (11) 

    # INPUT -   "equation"              - rovnica obsahujúca súradnice X,Y,Z zadefinovaná anonýmnou funkciou lambda, viď. DEFINE
    #           "X1,Y1...Z2"            - upravené súradnice, shifted endpoints (10)
    #           "*args"                 - dodatočné argumenty, pridať N (nmax) stupeň polynómu rho pri použití elemfun_n
    #           "homogenous"            - boolean, True - použijú sa elementárne funkcie pre homogénnnu prizmu, False - nehomogénna prizma

    # OUTPUT -  "triple_difference"     - vyčíslená hodnota zvolenej rovnice s použitím operátora trojitej diferencie 

    # DEFINE -  "equation"              - rovnica musí byť definovaná cez lambda (anonýmnu funkciu), kde parametrami funkcie sú súradnice XYZ a 
    #                                      tuple hodnoty vypočítaných element. funkcií pre danú prizmu v poradí v akom ich funkcia elemfun vracia 

    # ============================ F U N K C I A =================================== #

    if homogenous:                      # Predefinovanie názvu funkcie, použitie element. funkc. pre homogénnu či nehomogénnu prizmu 
        elemfun = elemfun_0
    else:
        elemfun = elemfun_n

    f1_element_val = elemfun(X2,Y2,Z2,*args)          # výpočet elementárnych funkcií pre rovnicu
    f1 = equation(X2, Y2, Z2,f1_element_val)    # výpočet rovnice na základe daných súradníc 

    f2_element_val = elemfun(X2,Y2,Z1,*args)
    f2 = equation(X2, Y2, Z1,f2_element_val)

    f3_element_val = elemfun(X2,Y1,Z2,*args)
    f3 = equation(X2, Y1, Z2,f3_element_val)

    f4_element_val = elemfun(X2,Y1,Z1,*args)
    f4 = equation(X2, Y1, Z1,f4_element_val)

    f5_element_val = elemfun(X1,Y2,Z2,*args)
    f5 = equation(X1, Y2, Z2,f5_element_val)

    f6_element_val = elemfun(X1,Y2,Z1,*args)
    f6 = equation(X1, Y2, Z1,f6_element_val)

    f7_element_val = elemfun(X1,Y1,Z2,*args)
    f7 = equation(X1, Y1, Z2,f7_element_val)

    f8_element_val = elemfun(X1,Y1,Z1,*args)
    f8 = equation(X1, Y1, Z1,f8_element_val)


    triple_difference = f1 - f2 - f3 + f4 - f5 + f6 + f7 - f8            # (11) -    trojitá diferencia funkcie (operátor)

    return triple_difference