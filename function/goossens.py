
def Goossens2Fukushima(prism,density):
    """
    Transforms density polynomial from Goossens conception (using z as depth coordinate) to Fukushima conception (using z as height coordinate).
    Currently only allows for lineary varying density. Thus, the function returns ValueError if given density of more (or less) than 2 coefficients. 

    Inputs:
    prism   ... list of 6 values,   the domain of the prism within the set coordinate system, x1,x2,y1,y2,z1,z2
    density ... list of 2 values,   the coefficients of the linear density polynomial, rho_0,rho_1

    Outputs:
    densityFukushima    ... list of 2 values, the coefficients of the linear density polynomial, rho_0,rho_1 transformed into using Fukushima's height coordinate system 
    """

    if len(prism) != 6: raise ValueError(f"Prism size incorrectly set, must be list of 6 values (the domain: x1,x2,y1,y2,z1,z2)")

    if len(density) != 2: raise ValueError(f"Density size set incorrectly, function currently only accepts linear variation of density, so list of coefficients: [rho0,rho1]")

    z_max = prism[5]
    rho_0G = density[0]
    rho_1G = density[1]

    # Simple transformation
    rho_0F = rho_0G + rho_1G * z_max
    rho_1F = - rho_1G

    return [rho_0F,rho_1F]