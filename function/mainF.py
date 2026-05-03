import numpy as np
import xarray as xr
import pygmt
import pandas as pd
from fukushima import Fukushima
from goossens import Goossens2Fukushima

def mainF(files, area, radius, mode, ref_sph=1727, density_mode="3d"):
    """
    Main function for computing gravitational effects using the Fukushima function along
    the goossens2fukushima to transform the density polynomial. WildPfeiffer formulas are used
    to transform input spherical coords into metric coords (tesseroids to right rectangular prisms)

    Inputs:
    files: list of DEM and density model
    area: chosen area to be evaluated
    radius: radius of the evaluation region
    mode: mode of operation (1: potential, 2: potential + acceleration, 3: potential + acceleration + tensor)
    ref_sph: reference sphere radius
    density_mode: mode for density evaluation ("constant", "2d", or "3d")

    Outputs:
    results: dictionary of calculated gravitational components
    stats: statistics from the evaluated region
    """
    # Reading input files
    DEM = xr.load_dataarray(files[0], engine='gmt', raster_kind='grid')
    RHO = xr.load_dataarray(files[1], engine='gmt', raster_kind='grid')
    GRAD = xr.load_dataarray(files[2], engine='gmt', raster_kind='grid')

    # Expand chosen area by one cell to compute the last row/col
    dlon = float(DEM.lon[1] - DEM.lon[0])
    dlat = float(DEM.lat[1] - DEM.lat[0])
    area = [area[0], area[1] + dlon, area[2], area[3] + dlat]

    # Convert radius from km to meters, and to an approximate cell radius for cutout
    radius_m = float(radius) * 1000.0
    ref_sph_m = float(ref_sph) * 1000.0
    phi0 = 0.5 * (area[2] + area[3])

    dy_cell = ref_sph_m * np.deg2rad(dlat)
    dx_cell = ref_sph_m * np.cos(np.deg2rad(phi0)) * np.deg2rad(dlon)
    cell_m = min(abs(dx_cell), abs(dy_cell))
    radius_cells = int(np.ceil(radius_m / cell_m))
    print(f"Radius in cells (approx): {radius_cells} cells for {radius_m/1000:.2f} km")

    # Cut out a region from DEM,RHO,GRAD of the chosen area + radius
    DEM_area, DEM_eval, RHO_area, RHO_eval, GRAD_area, GRAD_eval = cutout(DEM, RHO, GRAD, area, radius_cells)

    if density_mode not in ["constant", "2d", "3d"]:
        raise ValueError("density_mode must be 'constant', '2d', or '3d'")

    if density_mode == "constant":
        rho_mean = float(RHO_eval.mean())

    # Recalculate height in DEM - used to prevent negative heights
    height_offset = 1737.4 - ref_sph
    DEM_eval_m = (DEM_eval + height_offset) * 1000

    # Create Mercator projection of input files
    simfig(DEM_area + height_offset, cmap="viridis")
    simfig(RHO_area, cmap="plasma")
    simfig(GRAD_area, cmap="matlab/polar",ref_file=GRAD)

    # Recalculate spherical coordinates into local cartesian coordinates
    dx, dy, dz = tess2prism(DEM_eval_m, (ref_sph * 1000))

    # Pixel-registered loop setup
    lat_eval = DEM_eval_m.lat.values
    lon_eval = DEM_eval_m.lon.values
    lat_area = DEM_area.lat.values
    lon_area = DEM_area.lon.values

    # Output grid on cell centers (physics-correct for gridline registration)
    lat_c = (lat_area[:-1] + lat_area[1:]) / 2
    lon_c = (lon_area[:-1] + lon_area[1:]) / 2

    coords = {"lat": lat_c, "lon": lon_c}
    dims = ("lat", "lon")
    shape = (len(lat_c), len(lon_c))

    if mode >= 1:
        V = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
    if mode >= 2:
        gx = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        gy = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        gz = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
    if mode >= 3:
        gxx = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        gxy = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        gxz = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        gyy = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        gyz = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        gzz = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)

    x_center, y_center = cellcentre(dx, dy)

    lat_offset = np.where(np.isclose(lat_eval, lat_area[0]))[0][0]
    lon_offset = np.where(np.isclose(lon_eval, lon_area[0]))[0][0]

    total_eval = (len(lat_area) - 1) * (len(lon_area) - 1)
    max_prisms = (2 * radius_cells + 1) ** 2
    print(f"Max prisms per evaluation point (square window): {max_prisms}")
    print(f"Circular radius used: {radius_m/1000:.2f} km")

    # Loop over all cells in DEM_area (gridline registration: n-1 cells)
    for ia in range(len(lat_area) - 1):
        for ja in range(len(lon_area) - 1):

            ie = ia + lat_offset
            je = ja + lon_offset

            eval_index = ia * (len(lon_area) - 1) + ja + 1
            print(f"Evaluation point {eval_index}/{total_eval}")

            eval_point = [
                x_center[ie, je],
                y_center[ie, je],
                float(DEM_eval_m.isel(lat=ie, lon=je))
            ]

            V_sum = 0.0
            gx_sum = gy_sum = gz_sum = 0.0
            Gxx_sum = Gxy_sum = Gxz_sum = Gyy_sum = Gyz_sum = Gzz_sum = 0.0

            for i in range(ie - radius_cells, ie + radius_cells + 1):
                for j in range(je - radius_cells, je + radius_cells + 1):

                    if i < 0 or j < 0 or i >= len(lat_eval) - 1 or j >= len(lon_eval) - 1:
                        continue

                    cx = x_center[i, j]
                    cy = y_center[i, j]

                    # Circular radius check (in meters)
                    dxp = cx - eval_point[0]
                    dyp = cy - eval_point[1]
                    if (dxp * dxp + dyp * dyp) > (radius_m * radius_m):
                        continue

                    z_top = float(DEM_eval_m.isel(lat=i, lon=j))
                    z_bot = z_top - dz[i, j]

                    prism = [
                        cx - dx[i, j]/2, cx + dx[i, j]/2,
                        cy - dy[i, j]/2, cy + dy[i, j]/2,
                        z_bot, z_top
                    ]

                    if density_mode == "constant":
                        rho0 = rho_mean
                        rho1 = 0.0
                    elif density_mode == "2d":
                        rho0 = float(RHO_eval.isel(lat=i, lon=j))
                        rho1 = 0.0
                    else:
                        rho0 = float(RHO_eval.isel(lat=i, lon=j))
                        rho1 = float(GRAD_eval.isel(lat=i, lon=j))

                    density = Goossens2Fukushima(prism, [rho0, rho1])

                    if mode == 1:
                        V_val = Fukushima(prism, eval_point, density, mode=1)
                        V_sum += V_val
                    elif mode == 2:
                        V_val, g_val = Fukushima(prism, eval_point, density, mode=2)
                        V_sum += V_val
                        gx_sum += g_val[0]; gy_sum += g_val[1]; gz_sum += g_val[2]
                    elif mode == 3:
                        V_val, g_val, G_val = Fukushima(prism, eval_point, density, mode=3)
                        V_sum += V_val
                        gx_sum += g_val[0]; gy_sum += g_val[1]; gz_sum += g_val[2]
                        Gxx_sum += G_val[0]; Gxy_sum += G_val[1]; Gxz_sum += G_val[2]
                        Gyy_sum += G_val[3]; Gyz_sum += G_val[4]; Gzz_sum += G_val[5]

            if mode >= 2:
                gx_sum *= 1e5
                gy_sum *= 1e5
                gz_sum *= 1e5
            if mode >= 3:
                Gxx_sum *= 1e9
                Gxy_sum *= 1e9
                Gxz_sum *= 1e9
                Gyy_sum *= 1e9
                Gyz_sum *= 1e9
                Gzz_sum *= 1e9

            if mode >= 1:
                V[ia, ja] = V_sum
            if mode >= 2:
                gx[ia, ja] = gx_sum
                gy[ia, ja] = gy_sum
                gz[ia, ja] = gz_sum
            if mode >= 3:
                gxx[ia, ja] = Gxx_sum
                gxy[ia, ja] = Gxy_sum
                gxz[ia, ja] = Gxz_sum
                gyy[ia, ja] = Gyy_sum
                gyz[ia, ja] = Gyz_sum
                gzz[ia, ja] = Gzz_sum

    # Saving results in a dictionary for further use
    results = {}
    if mode >= 1:
        results["V"] = V
    if mode >= 2:
        results["gx"] = gx
        results["gy"] = gy
        results["gz"] = gz
    if mode >= 3:
        results["gxx"] = gxx
        results["gxy"] = gxy
        results["gxz"] = gxz
        results["gyy"] = gyy
        results["gyz"] = gyz
        results["gzz"] = gzz

    # Saving statistics in a list for further use
    stats = []

    # Also generate simple figures from the results
    if mode >= 1:
        stats.append(collect_stats("V [m^2/s^2]", V))
        simfig(V, cmap="plasma")

    if mode >= 2:
        stats.append(collect_stats("gx [mGal]", gx))
        stats.append(collect_stats("gy [mGal]", gy))
        stats.append(collect_stats("gz [mGal]", gz))
        simfig(gx, cmap="inferno")
        simfig(gy, cmap="inferno")
        simfig(gz, cmap="inferno")

    if mode >= 3:
        stats.append(collect_stats("Gxx [Eotves]", gxx))
        stats.append(collect_stats("Gxy [Eotves]", gxy))
        stats.append(collect_stats("Gxz [Eotves]", gxz))
        stats.append(collect_stats("Gyy [Eotves]", gyy))
        stats.append(collect_stats("Gyz [Eotves]", gyz))
        stats.append(collect_stats("Gzz [Eotves]", gzz))
        simfig(gxx, cmap="magma")
        simfig(gxy, cmap="magma")
        simfig(gxz, cmap="magma")
        simfig(gyy, cmap="magma")
        simfig(gyz, cmap="magma")
        simfig(gzz, cmap="magma")

    stats_df = pd.DataFrame(stats)
    print(stats_df.to_string(index=False))

    return results, stats


def simfig(map, cmap="viridis",ref_file=None, robust_pct=99):
    """
    Create a simple figure from a grid.

    Inputs:
    map: the grid to plot
    cmap: the colormap to use
    center_zero: whether to center the colormap at zero
    ref_grid: the reference grid for interpolation
    """

    if ref_file is None:
        fig = pygmt.Figure()
        with pygmt.config(MAP_FRAME_TYPE="fancy"):
            fig.grdimage(grid=map, projection="M15c", cmap=cmap)
            fig.basemap(frame=["a", "WSne"])
            fig.colorbar()
        fig.show()
    else:
        # Shorten colorbar using actual map values (robust symmetric range)
        vals = map.values
        vals = vals[np.isfinite(vals)]
        vmax = float(np.nanpercentile(np.abs(vals), robust_pct))
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1.0

        fig = pygmt.Figure()
        pygmt.makecpt(cmap=cmap, series=[-vmax, vmax])
        with pygmt.config(MAP_FRAME_TYPE="fancy", FORMAT_FLOAT_MAP="%.0f"):
            fig.grdimage(grid=map, projection="M15c", cmap=True)
            fig.basemap(frame=["a", "WSne"])
            fig.colorbar()
        fig.show()


def tess2prism(DEM, ref_sph_m):
    """
    Convert a spherical coordinate grid to a prism grid.
    Inputs:
    DEM: the digital elevation model
    ref_sph_m: the reference sphere radius

    Outputs:
    dx, dy, dz: metric dimensions of the prisms
    """
    lon = DEM.lon.values
    lat = DEM.lat.values

    dlon = float(lon[1] - lon[0])
    dlat = float(lat[1] - lat[0])

    dx = np.zeros(DEM.shape)
    dy = np.zeros(DEM.shape)
    dz = np.zeros(DEM.shape)

    for i in range(len(lat) - 1):
        for j in range(len(lon) - 1):
            lon1 = lon[j]
            lon2 = lon[j+1]
            lat1 = lat[i]
            lat2 = lat[i+1]

            r1 = ref_sph_m
            r2 = ref_sph_m + float(DEM.isel(lat=i, lon=j))

            dx[i, j], dy[i, j], dz[i, j] = WildPfeiffer(r1, r2, lon1, lon2, lat1, lat2)

    dx[-1, :] = dx[-2, :]
    dx[:, -1] = dx[:, -2]
    dy[-1, :] = dy[-2, :]
    dy[:, -1] = dy[:, -2]
    dz[-1, :] = dz[-2, :]
    dz[:, -1] = dz[:, -2]

    return dx, dy, dz


def cutout(DEM, RHO, GRAD, area, radius):
    """
    Extract a subset of the DEM within a specified area and radius.

    Inputs:
    DEM: the digital elevation model
    RHO: the density grid
    GRAD: the gravity gradient grid
    area: the bounding box coordinates (lon_min, lon_max, lat_min, lat_max)
    radius: the radius of the area to extract

    Outputs:
    DEM_area, DEM_eval, RHO_area, RHO_eval, GRAD_area, GRAD_eval
    """
    lon_min, lon_max, lat_min, lat_max = area

    dlon = float(DEM.lon[1] - DEM.lon[0])
    dlat = float(DEM.lat[1] - DEM.lat[0])

    ext_lon_min = lon_min - radius * dlon
    ext_lon_max = lon_max + radius * dlon
    ext_lat_min = lat_min - radius * dlat
    ext_lat_max = lat_max + radius * dlat

    DEM_eval = DEM.sel(lon=slice(ext_lon_min, ext_lon_max), lat=slice(ext_lat_min, ext_lat_max))
    DEM_area = DEM.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))

    RHO_eval = RHO.interp_like(DEM_eval, method='nearest')
    RHO_area = RHO.interp_like(DEM_area, method='nearest')
    GRAD_eval = GRAD.interp_like(DEM_eval, method='nearest')
    GRAD_area = GRAD.interp_like(DEM_area, method='nearest')

    return DEM_area, DEM_eval, RHO_area, RHO_eval, GRAD_area, GRAD_eval


def WildPfeiffer(r1, r2, lambda1, lambda2, phi1, phi2):
    """
    Calculate the metric dimensions of a prism based on its vertices.
    Inputs:
    r1, r2: radii of the two vertices
    lambda1, lambda2: longitudes of the two vertices
    phi1, phi2: latitudes of the two vertices

    Outputs:
    dx, dy, dz: metric dimensions of the prism
    """
    r0 = (r1 + r2) / 2
    phi0 = np.deg2rad((phi1 + phi2) / 2)
    dlambda = np.deg2rad(lambda2 - lambda1)
    dphi = np.deg2rad(phi2 - phi1)

    dx = r0 * dphi
    dy = r0 * np.cos(phi0) * dlambda
    dz = r2 - r1

    return dx, dy, dz


def cellcentre(dx, dy):
    """
    Calculate the center coordinates of each cell in a grid.
    Inputs:
    dx,dy: metric dimensions of prisms

    Outputs:
    centers of the evaluation points
    """
    x_center = np.zeros_like(dx)
    y_center = np.zeros_like(dy)

    for j in range(dx.shape[1]):
        for i in range(dx.shape[0]):
            if i == 0:
                x_center[i, j] = dx[i, j] / 2
            else:
                x_center[i, j] = x_center[i-1, j] + (dx[i-1, j]/2) + (dx[i, j]/2)

    for i in range(dy.shape[0]):
        for j in range(dy.shape[1]):
            if j == 0:
                y_center[i, j] = dy[i, j] / 2
            else:
                y_center[i, j] = y_center[i, j-1] + (dy[i, j-1]/2) + (dy[i, j]/2)

    return x_center, y_center


def collect_stats(label, grid):
    """
    Collect statistics from a grid.
    
    Inputs:
    label: name of the grid
    grid: grid variable

    Outputs:
    dictionary of stats of the grid
    """
    vals = grid.values
    vals = vals[np.isfinite(vals)]

    if vals.size == 0:
        return {
            "name": label,
            "min": np.nan,
            "max": np.nan,
            "mean": np.nan,
            "range": np.nan,
            "std": np.nan,
        }

    vmin = float(np.min(vals))
    vmax = float(np.max(vals))
    vmean = float(np.mean(vals))
    vstd = float(np.std(vals))
    vrange = vmax - vmin

    return {
        "name": label,
        "min": vmin,
        "max": vmax,
        "mean": vmean,
        "range": vrange,
        "std": vstd,
    }