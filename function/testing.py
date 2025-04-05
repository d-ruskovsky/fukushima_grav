import numpy as np
from fukushima import *
from tripledif import *

                # X      # Y    # Z
coor_prism = [-0.5,0.5, -1,1, -1.5,1.5] # !!!
coor_point = [1.5,3.0,4.5]
rho = [1]

# x1(0,0,0)
# x2()


V,g,G = Fukushima(coor_prism,coor_point,rho)

print("potencial =\n",V)
print("\nzrychlenie = \n",g)
print("\ntenzor =\n",G)

# VSTUPY 

# Newtonova gravitacna konstanta: 6.67430e-11 m**3 * kg*-1 * s**-2
# Rozmery prizmy (x, y, z): 1.0, 2.0, 3.0
# Hustota: 1.0 kg * m**-3
# Vypoctovy bod (xp, yp, zp): 1.5, 3.0, 4.5


# VYSTUPY 
# v     7.2025779027446533e-11 m**2 * s**-2 # potencial

# vx   -3.6702626407484486e-12 m * s**-2 # zrychlenia
# vy   -7.1486037410563400e-12 m * s**-2
# vz   -1.0322950752029843e-11 m * s**-2

# vxx  -1.8489526045069804e-12 s**-2 # tenzor
# vxy   1.1425556381054625e-12 s**-2
# vxz   1.6148227073209125e-12 s**-2

# vyy  -2.0256052043887284e-13 s**-2
# vyz   3.0899526561096312e-12 s**-2

# vzz   2.0515131249458494e-12 s**-2