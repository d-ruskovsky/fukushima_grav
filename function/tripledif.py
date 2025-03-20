import numpy as np

# Knižnica tripledif 
# Aktuálne obsahuje 2 funkcie - 
#                               evaluate_equation_0         - vyčíslenie rovnice (homogénna prizma - U_0)
#                               triple_dif                  - samotný operátor trojitej diferencie (Fukushima (11))

# TO DO: 
#   - doriešiť súradnice, kde sa hýbeme a jak, nerozumiem vzorcom (1),(6) a (10)
#   - vytvoriť funkcie na trojitú diferenciu aj pre nenulové funkcie U, "evaluate_equation_0" aktuálne používa parametre A,B..R použiteľné len pre homogénnu prizmu 


def evaluate_equation_0(rovnica,X,Y,Z):

    # funkcia pre evaluáciu zvolenej rovnice, aktuálne len pre rovnice homogénnej prizmy (vzorce Fukushima (23) až (25))

    # INPUT -   "rovnica"               - rovnica obsahujúca súradnice X,Y,Z (prípadne aj parametre elementárnych funkcií A,B,C,D,E,F a R) vo formáte string !! 
    #           "X,Y,Z"                 - trojrozmerné súradnice výpočtového bodu posunuté podľa (10) ??? nechápem, doriešiť 

    # OUTPUT -  "f"                     - vyčíslená hodnota zvolenej rovnice podľa zadaných súradníc 

    R = np.sqrt(X**2 + Y**2 + Z**2)                         # (25) -    euklidovská vzdialenosť

    D = np.log(X + R)                                       # (24) -    elementárne funkcie
    E = np.log(Y + R)
    F = np.log(Z + R)

    A = np.atan(Y*Z , X*R)
    B = np.atan(Z*X , Y*R)
    C = np.atan(X*Y , Z*R)

    f = eval(rovnica)

    return f



def triple_dif(rovnica,X1,Y1,Z1,X2,Y2,Z2):

    # funkcia trojitej diferencie aplikuje operátor podľa vzťahu Fukushima (11) na danú funkciu v stringu, zatiaľ platná len pre U_0 (homogénna prizma)!!!

    # INPUT -   "rovnica"               - rovnica obsahujúca súradnice X,Y,Z (prípadne aj parametre elementárnych funkcií A,B,C,D,E,F a R) vo formáte string !!
    #           "X1,Y1...Z2"            - trojrozmerné súradnice podľa článku, DORIEŠIŤ
    #           "evaluate_equation_0"   - funkcia vyčíslujúca vstupnú rovnicu so zadanými súradnicami a jej parametrami podľa Fukushima (24) a (25), čiže pre homogénnu prizmu U_0

    # OUTPUT -  "vysl"                  - vyčíslená hodnota zvolenej rovnice s použitím operátora trojitej diferencie 

    f1 = evaluate_equation_0(rovnica,X2, Y2, Z2)
    f2 = evaluate_equation_0(rovnica,X2, Y2, Z1)
    f3 = evaluate_equation_0(rovnica,X2, Y1, Z2)
    f4 = evaluate_equation_0(rovnica,X2, Y1, Z1)
    f5 = evaluate_equation_0(rovnica,X1, Y2, Z2)
    f6 = evaluate_equation_0(rovnica,X1, Y2, Z1)
    f7 = evaluate_equation_0(rovnica,X1, Y1, Z2)
    f8 = evaluate_equation_0(rovnica,X1, Y1, Z1)

    vysl = f1 - f2 - f3 + f4 - f5 + f6 + f7 - f8            # (11) -    trojitá diferencia funkcie

    return vysl