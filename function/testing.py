import numpy as np
from fukushima import *
from tripledif import *


coor_prism = [1,2,3,4,5,6]
coor_point = [20,30,40]
rho = [1]


V,g = Fukushima(coor_prism,coor_point,rho)

print("potencial =\n",V)
print("\nzrychlenie = \n",g)
