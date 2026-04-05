import logging
import numpy as np
import math

# Functions were made using these versions of above packages:
# numpy     2.4.2   - external, pyenv used
# logging   -       - Python standard library
# math      -       - Python standard library

# Create a logger for the whole library
logger = logging.getLogger("Fukushima")

def setup_logger(debug):
    """
    This function creates a logger for the below functions. Debugging is silenced unless set to True, in which case
    logging.info() is printed when running the main function. 
    """
    # Logger usage - info for debug mode, else only use warning
    logger.setLevel(logging.INFO if debug else logging.WARNING)
    # Formater for logging messeages
    formatter = logging.Formatter('%(asctime)s | %(funcName)s | %(message)s', datefmt='%H:%M:%S')
    # Prevent logger duplicity
    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)


def Fukushima(prism,point,density,mode=1,debug=False):
    """
    Function calculates the gravitational effect (potential V, acceleration g, tensor G) of a right rectangular 
    parallelepiped (prism) with a vertically varying density set by a arbitrary degree polynomial, based on an article 
    by Toshio Fukushima, DOI: 10.1093/gji/ggy317. 

    Inputs:
    prism   ... list of 6 values,   the domain of the prism within the set coordinate system, [x1,x2,y1,y2,z1,z2]
    point   ... list of 3 values,   the coordinates of the evaluation point [x,y,z]
    density ... list of n values,   the coefficients of the chosen density polynomial, rho_0,rho_1 ... rho_n
    mode    ... int,                1 (default) - compute gravitational potential only, 2 - gravitational potential + acceleration, 3 - potential + acceleration + tensor
    debug   ... boolean,            default set to False, True for using logger package to print each step of the calculation

    Outputs (based on mode set):
    V       ... float,              gravitational potential effect of the prism on the evaluation point
    g       ... list of float,      gravitational acceleration vector, [gx, gy, gz]
    G       ... list of float,      gravitational tensor, [Gxx, Gxy, Gxz, Gyy, Gyz, Gzz]
    """

    # Check if prism value set correctly
    if len(prism) != 6:
        raise ValueError(f"Prism must contain exactly 6 values [x1,x2,y1,y2,z1,z2], got {len(prism)}.")
    if not all(isinstance(v,(int,float)) for v in prism):
        raise ValueError(f"All prism values must be numeric")
    if prism[0] >= prism[1] or prism[2] >= prism[3] or prism[4] >= prism[5]:
        raise ValueError(f"Prism domain must satisfy x1<x2, y1<y2, z1<z2. Got x:[{prism[0]},{prism[1]}], y:[{prism[2]},{prism[3]}], z:[{prism[4]},{prism[5]}]")

    # Check if point value set correctly
    if len(point) != 3:
        raise ValueError(f"Point must contain exactly 3 values [x,y,z], got {len(prism)}.")
    if not all(isinstance(v, (int,float)) for v in point):
        raise ValueError(f"All point coordinates must be numeric")
    
    # Check if mode value set correctly
    if mode not in [1, 2, 3]:
        raise ValueError(f"Incorrect mode set: {mode}, allowed values: 1 (V), 2 (V+g), 3 (V+g+G). See documentation for further information.")
        
    # Set up logger for the main function
    setup_logger(debug)

    # Check if evaluation point is located within the prism - will not cause a crash, but may degrade results (see section 2.6)
    inside_x = prism[0] <= point[0] <= prism[1]
    inside_y = prism[2] <= point[1] <= prism[3]
    inside_z = prism[4] <= point[2] <= prism[5]
    if inside_x and inside_y and inside_z:
        logger.warning(f"Warning! Evaluation point within prism bounds.")
    
    # Calculated values according to set mode
    if mode==1:
        logger.info("Starting calculation of gravitational potential V")

    elif mode==2:
        logger.info("Starting calculation of gravitational potential V, acceleration g")
    
    elif mode==3:
        logger.info("Starting calculation of gravitational potential V, acceleration g, tensor G")
    
    # Shifting endpoints/vertices - (eq.10)
    # expectected prism domain in list - [x1,x2,y1,y2,z1,z2]
    #                                    [0, 1, 2, 3, 4, 5 ]
    # expected point coordinates in list - [x, y, z]
    #                                      [0, 1, 2]
    X1 = prism[0] - point[0]
    X2 = prism[1] - point[0]

    Y1 = prism[2] - point[1]
    Y2 = prism[3] - point[1]

    Z1 = prism[4] - point[2]
    Z2 = prism[5] - point[2]
    logger.info(f"Shifting endpoints - X1({X1}),X2({X2}),Y1({Y1}),Y2({Y2}),Z1({Z1}),Z2({Z2})")

    # Unpacking density polynomial degree
    # expected polynomial coefficients in 'density' list 
    N = len(density) - 1
    logger.info(f"Polynomial degree N set to {N}, polynomial coefficients: {density}")
    
    # Calculating gravitational effect of prism in evaluation point
    V = 0 # inicialize values for recursive computation

    if mode >= 2:
        # Gravitational acceleration components gx, gy, gz
        gx = 0
        gy = 0
        gz = 0

    if mode >= 3:
        # Gravitational force tensor components
        Gxx = 0
        Gxy = 0
        Gxz = 0
        Gyy = 0
        Gyz = 0
        Gzz = 0
    
    logger.info(f"Initial values set")

    for m in range(N + 1):
        logger.info(f"--Starting calculation, loop {m} of {N}")

        # Calculate the polynomial coefficient for current loop
        # c = cCoefficient(N, m, (prism[5]+Z2), density)
        c = cCoefficient(N, m, point[2], density)
        logger.info(f"Returned from function cCoefficient(), cm = {c[0]}, c'n = {c[1]}, c''n = {c[2]}, loop {m} of {N}")

        # Calculate the weight function for the current loop
        W = weightFunction(X1,X2,Y1,Y2,Z1,Z2,m,mode)
        logger.info(f"Returned from function weightFunction(), W = {W}, loop {m} of {N}")

        # Gravitational potential - (eq.13)
        V += c[0] * W["W"]
        logger.info(f"V = {V}, loop {m} of {N}")

        if mode >= 2:
            # Gravitational acceleration components (eq.16)
            gx += c[0] * W["Wx"]
            gy += c[0] * W["Wy"]
            gz += c[0] * W["Wz"] + c[1] * W["W"]
            logger.info(f"gx = {gx}, gy = {gy}, gz = {gz}")
        
        if mode >=3:
            # Gravitational force tensor components (eq.17)
            Gxx += c[0] * W["Wxx"]
            Gxy += c[0] * W["Wxy"]
            Gyy += c[0] * W["Wyy"]
            Gxz += c[0] * W["Wxz"] + c[1] * W["Wx"]
            Gyz += c[0] * W["Wyz"] + c[1] * W["Wy"]
            Gzz += c[0] * W["Wzz"] + 2 * c[1] * W["Wz"] + c[2] * W["W"]
            logger.info(f"Gxx={Gxx}, Gxy={Gxy}, Gxz={Gxz}, Gyy={Gyy}, Gyz={Gyz}, Gzz={Gzz}, loop {m} of {N}")

    # END of loop

    # END of main function
    if mode == 1:
        return V
    elif mode == 2:
        return V, [gx,gy,gz]
    elif mode == 3:
        return V, [gx,gy,gz], [Gxx,Gxy,Gxz,Gyy,Gyz,Gzz]

# Note: all functions below were made specifically to assist the Fukushima() function in evaluating the gravitational 
# effect of a right rectangular prism on an evaluation point. The functions therefore expect inputs to be already 
# accepted by the Fukushima() function, and can be prone to errors if used by themselves. 

def cCoefficient(N,m,z,density):
    """
    Calculates the value of polynomial coefficients and their 1st and 2nd derivatives for the specific loop set by Fukushima().
    Derivatives are only calculated if the polynomial is of a sufficiently high degree (1st derivative requires (N-m-1)>=0, 2nd derivative 
    requires (N-m-2)>=0). Furthermore the indexing of the derivatives is slightly changed to maintain consistency within the function.
    Indexes changed from article to function (article = function) -> m = j, n = m (only applies to derivatives). 

    Inputs:
    N       ... int, the degree of the density polynomial used
    m       ... int, currently calculated loop
    z       ... float, the height coordinate of the evaluation point
    density ... list of n values,   the coefficients of the chosen density polynomial, rho_0,rho_1 ... rho_n

    Outputs:
    [cm,cn_,cn__]      ... list, coefficient polynomial and its first and second derivatives with respect to z
    """

    cm = 0 # inicialize value of the coeficients for recursive computation
    cmj = 0

    cn_ = 0 # inicialize values of the 1st derivatives of polynomial coefficients
    cnm_ = 0

    cn__ = 0 # inicialize values of the 2nd derivatives of polynomial coefficients
    cnm__ = 0

    for j in range(N - m + 1):
        logger.info(f"Starting polynomial coefficient calculation, loop {j} out of (N - m) = {N - m} loops")
        
        # j-th coeficient as defined by (eq.15)
        cmj = math.comb((j + m), m) * (6.67430 * 10**-11) * density[(j+m)] # used G - Newtons constant of universal attraction, source B.Bucha Fyzikálna Geodézia 2023
        logger.info(f"Calculation of cmj = {cmj}, for loop {j} out of (N - m) {N - m}")

        # m-th coeficient as defined by (eq.14)
        cm += cmj * (z**j)
        logger.info(f"Calculation of cm = {cm}, for loop {j} out of (N - m) = {N - m}")

    # END of loop

    # Note: the article uses n-indexes for sumations of the coefficient derivatives. Here I am keeping m-indexing 
    # as in (eq.14,15) for consistency. 

    if (N - m - 1) < 0:
        logger.info(f"1st derivative of polynomial coefficient returning 0 as (N-n-1)={N-m-1}, must be at least 0")
        # Since the coefficient is already defined as 0, this line simply skips the calculation instead of redefining the variable
    else:
        for j in range(N - m - 1 + 1): # keeping +1 (Python) and -1 (eq.18) for context
            logger.info(f"Starting 1st derivative of the polynomial coefficient calculation, loop {j} out of (N - n - 1) = {N - m - 1}")

            # m-th coefficient 1st derivative (eq.19)
            cnm_ = (j + 1) * math.comb((m + j + 1),m) * (6.67430 * 10**-11) * density[(m+j+1)] # used G - Newtons constant of universal attraction, source B.Bucha Fyzikálna Geodézia 2023
            logger.info(f"Calculation of c'nm = {cnm_}, for loop {j} out of (N - m - 1) {N - m - 1}")

            # n-th coefficient 1st derivative (eq.18)
            cn_ += cnm_ * (z**j)
            logger.info(f"Calculation of c'n = {cn_}, for loop {j} out of (N - m - 1) = {N - m - 1}")

        # END of loop

    if (N - m - 2) < 0:
        logger.info(f"2nd derivative of polynomial coefficient returning 0 as (N-n-2)={N-m-2}, must be at least 0")
        # Since the coefficient is already defined as 0, this line simply skips the calculation instead of redefining the variable
    else:
        for j in range(N - m - 2 + 1): # keeping +1 (Python) and -2 (eq.18) for context
            logger.info(f"Starting 2nd derivative of the polynomial coefficient calculation, loop {j} out of (N - n - 2) = {N - m - 2}")

            # m-th coefficient 2nd derivative (eq.19)
            cnm__ = (j + 2) * (j + 1) * math.comb((m + j + 2),m) * (6.67430 * 10**-11) * density[(m+j+2)] # used G - Newtons constant of universal attraction, source B.Bucha Fyzikálna Geodézia 2023
            logger.info(f"Calculation of c''nm = {cnm__}, for loop {j} out of (N - m - 2) {N - m - 2}")

            # n-th coefficient 2nd derivative (eq.18)
            cn__ += cnm__ * (z**j)
            logger.info(f"Calculation of c'n = {cn__}, for loop {j} out of (N - m - 2) = {N - m - 2}")

        # END of loop

    return [cm,cn_,cn__]


def weightFunction(X1,X2,Y1,Y2,Z1,Z2,n,mode):
    """
    Works as a bridge between several functions, picking the correct potential function and 
    returning the function after the triple difference operator has been applied to it.

    Inputs:
    X1,X2,Y1,Y2,Z1,Z2   ... shifted endpoints, the domain of the prism shifted by the coordinates of the evaluation point (eq.10)
    n                   ... currently calculated loop
    mode                ... int, 1 (default) - compute gravitational potential only, 2 - gravitational potential + acceleration, 3 - potential + acceleration + tensor

    Outputs:
    result              ... dictionary of weight functions, potential function W is always calculated (mode 1), the result 
                            also includes Wx,Wy,Wz (weight functions of acceleration vector g) for mode 2, and Wxx,Wxy,Wxz,Wyy,Wyz,Wzz 
                            (weight functions of the force tensor G) for mode 3
    """

    if n == 0:
        logger.info(f"Used potential functions for homogenous prism, n = {n}")
        # Potential functions for a homogenous prism (eq.23)

        # Gravitational potential V
        def U(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return -((X**2 * e["A"]) + (Y**2 * e["B"])\
                    + (Z**2 * e["C"]))/2\
                    + Y * Z * e["D"][0] + Z * X * e["E"][0] + X * Y * e["F"]
        
        # Gravitational acceleration vector g
        def Ux(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return - X * e["A"] + Y * e["F"] + Z * e["E"][0]
        
        def Uy(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return X * e["F"] - Y * e["B"] + Z * e["D"][0]
        
        def Uz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return X * e["E"][0] + Y * e["D"][0] - Z * e["C"]
        
        # Gravitational force tensor G
        def Uxx(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return -e["A"]
        
        def Uxy(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return e["F"]
        
        def Uxz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return e["E"][0]
        
        def Uyy(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return -e["B"]
        
        def Uyz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return e["D"][0]
        
        def Uzz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return -e["C"]
    
    if n >= 1:
        logger.info(f"Used potential functions for nonhomogenous prismn, n = {n}")
        # Potential function for a non-homogenous prism

        # Gravitational potential (eq.26)
        def U(X,Y,Z): 
            e = elementaryFunction(X, Y, Z, n)
            return -((Z**(n+2) * e["C"])/(n + 2))\
                     + ((Z**(n+1) * (Y * e["D"][1] + X * e["E"][1]))/(n + 1))\
                     - ((Y * e["D"][n+2] + X * e["E"][n+2])/((n + 1) * (n + 2)))
        
        # Gravitational acceleration vector (eq.33)
        def Ux(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return ((Z**(n + 1) * e["E"][0] - e["E"][n+2])/(n + 1))
        
        def Uy(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return ((Z**(n + 1) * e["D"][0] - e["D"][n+2])/(n + 1))
        
        def Uz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return -Z**(n + 1) * e["C"] + Z**n * (Y * e["D"][0] + X * e["E"][0])
        
        # Gravitational force tensor (eq.33)
        def Uxx(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return X * e["E"][n]
        
        def Uxy(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return e["R"][n]
        
        def Uxz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return Z**n * e["E"][0]
        
        def Uyy(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return Y * e["D"][n]
        
        def Uyz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return Z**n * e["D"][0]
        
        def Uzz(X,Y,Z):
            e = elementaryFunction(X, Y, Z, n)
            return -(n + 1) * Z**n * e["C"] + n * Z**(n - 1) * (Y * e["D"][0] + X * e["E"][0])
        
    # Weight function expressed and rewritten using triple difference (eq.9)
    result = {"W": tripleDifference(U,X1,X2,Y1,Y2,Z1,Z2)}
    logger.info(f"Weight function W for potential evaluated: {result['W']}")

    if mode >= 2:
        # 1st derivatives of the weight function (acceleration vector g) (eq.20)
        result["Wx"] = -tripleDifference(Ux,X1,X2,Y1,Y2,Z1,Z2)
        result["Wy"] = -tripleDifference(Uy,X1,X2,Y1,Y2,Z1,Z2)
        result["Wz"] = -tripleDifference(Uz,X1,X2,Y1,Y2,Z1,Z2)
        logger.info(f"Weight functions Wx,Wy,Wz for acceleration vector evaluated: {result['Wx']},{result['Wy']},{result['Wz']}")

    if mode >= 3:
        # 2nd derivatives of the weight function (force tensor G) (eq.20)
        result["Wxx"] = tripleDifference(Uxx,X1,X2,Y1,Y2,Z1,Z2)
        result["Wxy"] = tripleDifference(Uxy,X1,X2,Y1,Y2,Z1,Z2)
        result["Wxz"] = tripleDifference(Uxz,X1,X2,Y1,Y2,Z1,Z2)
        result["Wyy"] = tripleDifference(Uyy,X1,X2,Y1,Y2,Z1,Z2)
        result["Wyz"] = tripleDifference(Uyz,X1,X2,Y1,Y2,Z1,Z2)
        result["Wzz"] = tripleDifference(Uzz,X1,X2,Y1,Y2,Z1,Z2)
        logger.info(f"Weight functions Wxx,Wxy,Wxz,Wyy,Wyz,Wzz for force tensor evaluated: {result['Wxx']},{result['Wxy']},{result['Wxz']},{result['Wyy']},{result['Wyz']},{result['Wzz']}")

    return result


def elementaryFunction(X,Y,Z,n):
    """
    Calculates the elementary functions of the specific coordinate combination for a set degree(loop) n. 
    Elementary functions A,B,C,F are evaluated as floats, given they are not computed recursively. 
    Elementary functions D,E,R are evaluated as lists of floats, starting at [0] for n=0 (homogenous prism) 
    and ending at [n+2] as required by potential functions. These can then be called easily, for example
    elementary["D"][2] will return the value of the elementary function D for n=2. 

    Inputs:
    X,Y,Z       ... float, coordinates
    n           ... current polynomial degree index

    Outputs:
    elementary  ... dict, comprised of floats (A,B,C,F) and lists (D,E,R) evaluated for the specific coordinates and degree/loop (n+2)

    """
    logger.info(f"Calculating elementary functions for X = {X}, Y = {Y}, Z = {Z}. Loop {n}")

    # Squared sum of the horizontal coordinates (eq.31)
    S = X**2 + Y**2

    # Initial values for elementary functions (eq.36)
    A = atan3(X,Y,Z)
    B = atan3(Y,Z,X)
    C = atan3(Z,X,Y)
    F = logsum(Z,X,Y)
    D = [None] * (n + 3) # Note: creating an empty list to fill using recursive computation
    E = [None] * (n + 3)
    R = [None] * (n + 3)

    D[0] = logsum(X,Y,Z)
    E[0] = logsum(Y,Z,X)
    R[0] = np.sqrt(X**2 + Y**2 + Z**2 + (10**-150)) # Note: value 10**-150 is defined as omega (eq.39) - using the value directly

    logger.info(f"Initial values of elementary functions: A:{A},B:{B},C:{C},D:{D[0]},E:{E[0]},F:{F},R:{R[0]}")

    # Recursive computation uses elementary functions of a homogenous prism as initial values (eq.28,30)
    D[1] = D[0]
    E[1] = E[0]
    R[1] = R[0]

    # Second degree elementary functions also have specific equations (eq.28,30)
    D[2] = Y*B - X*F
    E[2] = X*A - Y*F
    R[2] = (Z*R[0] - S*F)/2

    logger.info(f"Elementary functions for n=2, D[2]:{D[2]},E[2]:{E[2]},R[2]:{R[2]}")

    logger.info(f"Starting recursive computation of elementary functions")
    for m in range(3,n + 2 + 1): # + 2 for the need of _n+2, + 1 for Python - check later if correct

        # Note: this loop starts at 3, given values for n=2 are already evaluated. The loop therefore 
        # starts at 3 and continues to calculate recursively until n+2 as required by potential
        # functions (eq.26) and (eq.33)

        R[m] = (Z**(m-1) * R[0] - (m - 1) * S * R[m-2])/m   # (eq.29)
        D[m] = -Y**2 * D[m-2] - X * R[m-2]                  # (eq.27)
        E[m] = -X**2 * E[m-2] - Y * R[m-2]                  # (eq.27)

        logger.info(f"Elementary functions for n={m}, D[{m}]:{D[m]},E[{m}]:{E[m]},R[{m}]:{R[m]}")

    # END of loop
        
    return {"R":R, "A":A, "B":B, "C":C, "D":D, "E":E, "F":F}


def tripleDifference(function,X1,X2,Y1,Y2,Z1,Z2):
    """
    Evaluates a potential function set by weightFunction() for a specific vertex set by
    the vertexes coordinates (shifted endpoints), then applies the triple difference (eq.11) operator to it.

    Inputs:
    function            ... object, the equation/potential function, set by weightFunction()
    X1,X2,Y1,Y2,Z1,Z2   ... float, shifted endpoints, the domain of the prism shifted by the coordinates of the evaluation point (eq.10)

    Outputs:
    triple difference   ... float, potential function's alternating sum over all 8 vertices
    """
    
    logger.info(f"Triple difference applied for function {function},X1({X1}),X2({X2}),Y1({Y1}),Y2({Y2}),Z1({Z1}),Z2({Z2}")

    return function(X2,Y2,Z2) - function(X2,Y2,Z1) - function(X2,Y1,Z2) + function(X2,Y1,Z1) \
            - function(X1,Y2,Z2) + function(X1,Y2,Z1) + function(X1,Y1,Z2) - function(X1,Y1,Z1)

# Note: The author mentions in section 2.6, that evaluation near or inside verices or surface of the prism may lead to 
# degradation of the results. The bellow functions will prevent a crash in calculation, however the resulting values 
# might be damaged. 

def atan3(X,Y,Z):
    """
    Calculates the arctangent of specified X,Y,Z coordinates, while 
    using the branch structure for returning a 0 if the xi coordinate is also 0 (eq.37)

    Inputs:
    X,Y,Z   ... float, coordinates

    Outputs:
    atan3   ... float, the arctangent value
    """
    # Possibly worth noting that the coordinates are written as xi, eta and zeta in the article.
    # Or maybe not worth noting, but note it I shall nonetheles. 

    logger.info(f"Evaluating the arcus tangents of X:{X},Y:{Y},Z:{Z}")

    if X == 0:
        return 0
    else:
        return np.arctan((Y * Z) / (X * np.sqrt(X**2 + Y**2 + Z**2 )))  
    
def logsum(X,Y,Z):
    """
    Calculates the natural logarithm of specified X,Y,Z coordinates, while 
    using the branch structure that selects among equivalent forms to prevent numerical errors. (eq.38)

    Inputs:
    X,Y,Z       ... float, coordinates

    Outputs:
    logsum      ... float, the natural logarithm value
    """
    # Value 10**-150 in the article is defined as omega "tiny positive constant". Decided to use the value 
    # in the function directly instead of creating a separate variable to prevent needless re-defining 
    # of the variable at each loop. 

    logger.info(f"Evaluating the natural logarithm of X:{X},Y:{Y},Z:{Z}")

    if X > 0:
        return np.log(X + np.sqrt(X**2 + Y**2 + Z**2))
    elif X == 0:
        return (1/2) * np.log(Y**2 + Z**2 + (10**-150))
    else:
        return np.log((Y**2 + Z**2 + (10**-150)) / (np.sqrt(X**2 + Y**2 + Z**2) - X))