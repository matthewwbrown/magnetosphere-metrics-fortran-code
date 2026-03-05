import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os

def compute_ISA_timeseries(folder, region, file_range, dt_minutes=5):
    (x1, x2, y1, y2, z1, z2) = region 
    results = []  #Define empty array for the ISA at every timeframe
    times = []    #Define empty array for the timeframe that is being inspected
    x_width = abs(x2 - x1) 
    y_width = abs(y2 - y1) 
    z_width = abs(z2 - z1)   

    for i in range(file_range[0], file_range[1] + 1): #Iterates through file indices
        filename = os.path.join(folder, f"rate_{i}.nc") 
        if not os.path.exists(filename):
            continue  #If a filename is missing, then skip over

        R = xr.open_dataset(filename, engine="netcdf4") #loads netCDF file using xarray

        #mask a spherical region of 2 Earth radii around the Earth
        R = R.where(np.sqrt(R.x**2 + R.y**2 + R.z**2) > 2)

        rate = R["rate"] #grabs the 3D array named rate

        # Extract subvolume defined by region
        rate_sel = rate.sel(
            x=slice(x1, x2),
            y=slice(y1, y2),
            z=slice(z1, z2)
        )

        # Define the grid spacing (uniform)
        dx = float(rate.x[1] - rate.x[0])
        dy = float(rate.y[1] - rate.y[0])
        dz = float(rate.z[1] - rate.z[0])

        integral = np.nansum(rate_sel.values) * dx * dy * dz #Riemann sum over the rate with respect to volume, ignoring any NaN values
        mean_slip = integral / (x_width*y_width*z_width) #Normalise to the volume of the defined region to get an averge across the whole region

        results.append(mean_slip)
        times.append(dt_minutes * (i - file_range[0])) #Adjust for first file corresponds to file 0 

    return np.array(times), np.array(results)


def compute_unsigned_flux_timeseries(folder, x_plane, file_range, dt_minutes=5):
    flux_vals = []
    times = []

    for i in range(file_range[0], file_range[1] + 1): #Iterates through file indices
        filename = os.path.join(folder, f"output_{i}.nc")
        if not os.path.exists(filename):
            continue #If a filename is missing, then skip over

        ds = xr.open_dataset(filename, engine="netcdf4")

        #mask a sphericval region of 2 Earth radii around the Earth
        ds = ds.where(np.sqrt(ds.x**2 + ds.y**2 + ds.z**2) > 2)

        # Extract Bx at the x location closest to x_plane
        Bx_plane = ds["bx"].sel(x=x_plane, method="nearest")

        # Define the grid spacing (uniform)
        dy = float(ds.y[1] - ds.y[0])
        dz = float(ds.z[1] - ds.z[0])

        # Unsigned flux per plane
        flux = np.nansum(np.abs(Bx_plane.values)) * dy * dz

        flux_vals.append(flux)
        times.append(dt_minutes * (i - file_range[0])) #Adjust for first file corresponds to file 0 

        ds.close()

    return np.array(times), np.array(flux_vals)


