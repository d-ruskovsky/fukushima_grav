import rasterio
from rasterio.windows import Window

# Súbor na stríhanie (clipping) rastrov .tif podľa pixelov


# Vstupné a výstupné súbory =============================
INPUT_FILE = 'grid/extracted_100m.tif'

OUTPUT_FILE = 'grid/clipped/100m.tif'

# Clipping parametre ====================================
x_offset = 450 # poč. stlpec
y_offset = 450 # poč. riadok
width = 50 # rozmer clipped rastra (v pixeloch)
height = 50 # rozmer clipped rastra (v pixeloch)

# Clipping ==============================================
with rasterio.open(INPUT_FILE) as src:
    # Definovať okno ktorým sa raster clipne
    window = Window(x_offset,y_offset,width,height)

    # Vytvorenie profilu pre nový raster
    profile = src.profile
    profile.update({
        "height": height,
        "width": width,
        "transform": src.window_transform(window)
    })

    # Uložiť clipnutý raster 
    with rasterio.open(OUTPUT_FILE, "w", **profile) as dst:
        dst.write(src.read(window = window))
