import rasterio 

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Najvhodnejšie farby - 'viridis','magma','inferno','plasma','cividis'

def vykreslenie_tif(INPUT_FILE,colorbar_label,farby):
# Vstupné dáta ==================================
    INPUT_FILE = rasterio.open(INPUT_FILE)

    # Prečítanie výšok rastra
    raster_data = INPUT_FILE.read(1)
    # Extrakcia rozmerov, súradnice S-JTSK
    bounds = INPUT_FILE.bounds

    # Matplotplib.pyplot premenné pre grafiku
    fig, ax = plt.subplots(figsize=(10,10))
    cax = ax.imshow(raster_data, extent=(bounds.left, bounds.right, bounds.bottom, bounds.top), cmap = farby)
    # Definovanie mimorámových údajov 
    ax.set_title("Záujmové územie (Stredné Slovensko)",fontsize=16)

    ax.set_xlabel("Y [m]",fontsize = 12)
    ax.set_ylabel("X [m]",fontsize = 12)

    # Zrušenie vedeckej notácie 
    ax.ticklabel_format(style="plain",axis="both")

    # Pootočenie súradníc osi y (JTSK-X)
    ax.tick_params(axis='y', labelrotation = 90)

    # Colorbar
    colorbar = fig.colorbar(cax, ax=ax, shrink=0.7, fraction=0.05, pad=0.02)
    colorbar.set_label(colorbar_label , fontsize = 12)

    # Nastavenie gridu
    plt.grid(visible=True)

    # Dodatočný text
    fig.text(0.98, 0.05, "Súradnicový systém: S-JTSK (JTSK03)", fontsize=12, color="black",ha="right")

    # Zobrazenie plotu
    plt.tight_layout() # Krajší layout
    plt.show()