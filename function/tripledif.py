import numpy as np

# Knižnica tripledif 

# Aktuálne obsahuje 3 funkcie - 
#                               elemfun_0              - elementárne funkcie (homogénna prizma)
#                               elemfun_n              - elementárne funkcie (nehomogénna prizma)
#                               triple_dif             - operátor trojitej diferencie (Fukushima (11))

# TO DO: 
#            urobiť funkciu pre elemfun_n 

def elemfun_0(X,Y,Z):

    # funkcia pre výpočet elementárnych funkcií A,B,C,D,E,F (24) v prípade homogénnej prizmy 

    # INPUT -   "X,Y,Z"                 - trojrozmerné súradnice bodu posunuté podľa (10) 

    # OUTPUT -  "A,B,C,D,E,F"           - elementárne funkcie podľa (24) v prípade homogénnej prizmy  

    R = np.sqrt(X**2 + Y**2 + Z**2)                         # (25) -    euklidovská vzdialenosť

    D = np.log(X + R)                                       # (24) -    elementárne funkcie
    E = np.log(Y + R)
    F = np.log(Z + R)
    
    A = np.arctan((Y*Z) / (X*R))
    B = np.arctan((Z*X) / (Y*R))
    C = np.arctan((X*Y) / (Z*R))

    return A,B,C,D,E,F


def elemfun_n(X,Y,Z):

    return 0



def triple_dif(equation,X1,Y1,Z1,X2,Y2,Z2,homogenous=True):

    # funkcia trojitej diferencie aplikuje operátor podľa vzťahu Fukushima (11) 

    # INPUT -   "equation"              - rovnica obsahujúca súradnice X,Y,Z zadefinovaná anonýmnou funkciou lambda, viď. DEFINE
    #           "X1,Y1...Z2"            - upravené súradnice, shifted endpoints (10)
    #           "homogenous"            - boolean, True - použijú sa elementárne funkcie pre homogénnnu prizmu, False - nehomogénna prizma

    # OUTPUT -  "vysl"                  - vyčíslená hodnota zvolenej rovnice s použitím operátora trojitej diferencie 

    # DEFINE -  "equation"              - rovnica musí byť definovaná cez lambda (anonýmnu funkciu), kde parametrami funkcie sú súradnice XYZ a 
    #                                      tuple hodnoty vypočítaných element. funkcií pre danú prizmu v poradí v akom ich funkcia elemfun vracia 

    if homogenous:                      # Predefinovanie názvu funkcie, použitie element. funkc. pre homogénnu či nehomogénnu prizmu 
        elemfun = elemfun_0
    else:
        elemfun = elemfun_n

    f1_element_val = elemfun(X2,Y2,Z2)          # výpočet elementárnych funkcií pre funkciu
    f1 = equation(X2, Y2, Z2,f1_element_val)    # výpočet funkcie na základe daných súradníc 

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


    vysl = f1 - f2 - f3 + f4 - f5 + f6 + f7 - f8            # (11) -    trojitá diferencia funkcie

    return vysl