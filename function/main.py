import time
import numpy as np
import rasterio
from fukushima import *
from vykreslenie import *

start_time_raster = time.time()  # Začiatok odpočtu času

# ================== Vytváranie rastrov ================== 
input_raster = "grid/clipped/100m.tif" # zvolenie vstupného rastra
output_rasters = [                      # zvolenie výsledných rastrov
    "grid/final/V_0_100m.tif",
    "grid/final/gx_0_100m.tif",
    "grid/final/gy_0_100m.tif",
    "grid/final/gz_0_100m.tif"
]

created_rasters = [] # Prázdny array pre názvy vytvorených rastrov
with rasterio.open(input_raster) as src:
    # Prevzatie metadát z pôvodného rastra
    profile = src.profile

    # Načítanie S-JTSK 
    bounds = src.bounds
    top_left = (abs(bounds.left),abs(bounds.top))

    # Prečítanie výšok pixelov 
    h = src.read(1)

    # Nastavenie nulových hodnôt pre dáta v novom rastri
    profile.update({
        "dtype": "float32",  # zmeň datatyp ak potrebné
        "nodata": 0         # no data
    })

    # Vytvorenie prázdneho array rovnakej veľkosti ako pôvodný raster
    empty_data = np.full((profile['count'], profile['height'], profile['width']),
                          profile['nodata'], dtype=profile['dtype'])

    # Napíš nový raster
    for output_raster in output_rasters:
        with rasterio.open(output_raster, "w", **profile) as dst:
            dst.write(empty_data)

        # Ukladanie názvov rastrov pre ďalšie použitie (vykreslenie)
        created_rasters.append(output_raster)

# Ukončenie vytvárania rastrov
print(f"Vytvorenie rastrov trvalo: {time.time() - start_time_raster:.2f} sekúnd")

# =================== Výpočet gravitačného účinku ================

start_time_fukushima = time.time() # začiatok odpočtu času pre výpočet fukushima

rho = [1] # polynóm hustoty 

rows, columns = src.shape

out_V = np.full((rows, columns), profile['nodata'], dtype=np.float32)
out_gx = np.full((rows, columns), profile['nodata'], dtype=np.float32)
out_gy = np.full((rows, columns), profile['nodata'], dtype=np.float32)
out_gz = np.full((rows, columns), profile['nodata'], dtype=np.float32)

V0 = 0
g0 = np.array([0,0,0])

for k in range(rows):
    for l in range(columns):

        V0 = 0
        g0 = np.array([0,0,0])
        
        # coor_point = [(top_left[1] + 50) + (k*100),(top_left[0] + 50) + (l*100),h[k,l]]
        coor_point = [50 + (k*100), 50 + (l*100),h[k,l]]
        for i in range(rows):
            for j in range(columns):
            
                # coor_prism = [top_left[1] - i*100,top_left[1] - (i+1)*100, \
                #              top_left[0] - j*100, top_left[0] - (j+1)*100,\
                #                  0,h[i,j]]
                
                coor_prism = [i*100,(i+1)*100, \
                              j*100,(j+1)*100,\
                                  0,h[i,j]]
                V,g,_ = Fukushima(coor_prism,coor_point,rho)
                V0 = V0 + V 
                g0 = g0 + g

        out_V[k,l] = V0
        out_gx[k,l] = g0[0]
        out_gy[k,l] = g0[1]
        out_gz[k,l] = g0[2]
        print(f"Výpočet pixelu [{k},{l}] trval: {time.time() - start_time_fukushima:.2f} sekúnd")

        

# Option 1: Write the out_V array to the V_0_100m.tif raster file
with rasterio.open("grid/final/V_0_100m.tif", "r+") as dst:
    # Write the computed results to band 1 of the output raster.
    dst.write(out_V, 1)
with rasterio.open("grid/final/gx_0_100m.tif", "r+") as dst:
    dst.write(out_gx, 1)
with rasterio.open("grid/final/gy_0_100m.tif", "r+") as dst:
    dst.write(out_gy, 1)
with rasterio.open("grid/final/gz_0_100m.tif", "r+") as dst:
    dst.write(out_gz, 1)

# Ukončenie výpočtu rastrov
print(f"Výpočet gravitačného účinku trval: {time.time() - start_time_fukushima:.2f} sekúnd")

vykreslenie_tif(input_raster,"Nadmorská výška [m]",'viridis')
vykreslenie_tif(created_rasters[0],"Gravitačný potenciál [m^2 s^-2]",'cividis')
vykreslenie_tif(created_rasters[1],"Gravitačné zrýchlenie (x-zložka) [m s^-2]",'plasma')
vykreslenie_tif(created_rasters[2],"Gravitačné zrýchlenie (y-zložka) [m s^-2]",'plasma')
vykreslenie_tif(created_rasters[3],"Gravitačné zrýchlenie (z-zložka) [m s^-2]",'plasma')