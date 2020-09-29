import xarray as xr
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import cmocean

from Tests2 import calc_trajectory
from Tests2 import vertical_prof_trajectory
from Tests2 import plot_map_background
from Tests2 import XYT_NonDepth_1VAR_trajectory
from Tests2 import Depth_2VAR_trajectory
from Tests2 import XYU_NonDepth_2VAR_trajectory

###############################################################################

# Cross (Butterfly) Sectioning (Saildrone 1)
trajectory = {'datetime': [dt.datetime(2003, 9, 5, 0, 0, 0)],
              'lon': [-105], 'lat': [15.5]}
##                          spd  hours direction   attach to this trajectory
trajectory = calc_trajectory(1.25, 144, -90.0, prev_trajectory=trajectory)
trajectory = calc_trajectory(1.25, 144, 90.0, prev_trajectory=trajectory)

# (Saildrone 2)
trajectory2 = {'datetime': [dt.datetime(2003, 9, 5, 0, 0, 0)],
              'lon': [-107.825], 'lat': [12.5]}
##                          spd  hours direction   attach to this trajectory
trajectory2 = calc_trajectory(1.25, 144, 0.0, prev_trajectory=trajectory2)
trajectory2 = calc_trajectory(1.25, 144, 180.0, prev_trajectory=trajectory2)

###############################################################################

# Import data

# Open dataset via xarray
ds = xr.open_mfdataset('D:/CFSv2/Tests/*.nc')

# Convert time from object to datetime64
time = ds.variables['time'][:]
datetimeindex = ds.indexes['time'].to_datetimeindex()  # converts to python friendly 'datetime'
ds['time'] = datetimeindex

# Other Slices
start_time = '2003-09-05 00:00:00'
end_time = '2003-09-17 23:00:00'
time_slice = slice(start_time, end_time)
lat_slice = slice(6, 20)
lon_slice = slice(-115, -95)
depth_slice = slice(0, 150)
depth = 0

# Times in hours (Product 288hrs/24hrs = 12days (based on time slices)
Vtime = np.arange(0, 289, 1)

########################################################################

# Plot 1 (Temp and Salinity)

# Read in Vertical Profiles (var, data, ...)
F = vertical_prof_trajectory('temp', ds, depth_slice, lon_slice, lat_slice, time_slice, trajectory)
G = vertical_prof_trajectory('salt', ds, depth_slice, lon_slice, lat_slice, time_slice, trajectory)
H = vertical_prof_trajectory('temp', ds, depth_slice, lon_slice, lat_slice, time_slice, trajectory2)
I = vertical_prof_trajectory('salt', ds, depth_slice, lon_slice, lat_slice, time_slice, trajectory2)

plt.style.use('dark_background')

# Day 1 Data
paramt = 'temp'
params = 'salt'
depth2 = 0.5

# Start and end time for day 1
ST = '2003-09-05 00:00:00'
ET = '2003-09-05 23:00:00'
time_slice1 = slice(ST, ET)

# Get data, selecting time, depth, lat/lon slice
SST = ds[paramt].sel(st_ocean=depth2, time=time_slice1)
SSS = ds[params].sel(st_ocean=depth2, time=time_slice1)

# Specify longitude and latitude values for chosen domain
lons1 = SST.xt_ocean.values
lats1 = SST.yt_ocean.values

# average the data by day
avg_SST = SST.sum(dim='time') / 24  # 24 for one day
avg_SSS = SSS.sum(dim='time') / 24

################################################################################

# Set Up Figure
plot_area = [-115, -95, 6, 20]  # [lon1, lon2, lat1, lat2], lon is '-280'-'90' for CFSv2 model.

# levels and cbar_ticks sets up color scales
cmap00 = plt.cm.jet
cmap0 = LinearSegmentedColormap.from_list('custom', cmap00(np.linspace(0.10, 1.0, 28)))
cmap = cmap0
levels = MaxNLocator(nbins=28).tick_values(10.0, 34.0)
norm = BoundaryNorm(levels, cmap0.N, clip=True)

fig = plt.figure(figsize=(8.5, 7.5))

#################################################################################

# SST Map
ax1 = fig.add_subplot(3, 2, 1)
map1 = plot_map_background(plot_area=plot_area, ax=ax1)
H1 = plt.pcolormesh(lons1, lats1, avg_SST, cmap=cmap, norm=norm)

# Plots tracks on top graph
plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')  # Change if distance changes

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H1)
cbar.set_ticks(np.arange(12, 34, 4))

ax1.set_title('Temperature')

# SST Transect SD1
ax3 = fig.add_subplot(3, 2, 3)

# See functions from Tests2
H3 = plt.contourf(Vtime, F['Depths'], F['Verts'], levels=levels, cmap=cmap, norm=norm, extend='both')
plt.gca().invert_yaxis()  # Invert cause depths

cbar = plt.colorbar(H3)
cbar.set_ticks(np.arange(12, 34, 4))

ax3.set_ylim([100, 0])
ax3.set_xlabel('Time [h]')
ax3.set_ylabel('Depth [m]')
ax3.set_title('Temperature [A]')

#################################################################################

# SSS Map
cmap00 = cmocean.cm.haline
cmap0 = LinearSegmentedColormap.from_list('custom', cmap00(np.linspace(0.0, 1.0, 15)))
cmap_sss = cmap0
levels = MaxNLocator(nbins=15).tick_values(32.5, 36.6)
norm_sss = BoundaryNorm(levels, cmap0.N, clip=True)

ax2 = fig.add_subplot(3, 2, 2)
map2 = plot_map_background(plot_area=plot_area, ax=ax2)
H2 = plt.pcolormesh(lons1, lats1, avg_SSS, cmap=cmap_sss, norm=norm_sss)

plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H2)
cbar.set_ticks(np.arange(33, 37, 1))

ax2.set_title('Salinity')

## SSS Transect SD 1
ax4 = fig.add_subplot(3, 2, 4)
H4 = plt.contourf(Vtime, G['Depths'], G['Verts'], levels=levels, cmap=cmap_sss, norm=norm_sss, extend='both')
plt.gca().invert_yaxis()

cbar = plt.colorbar(H4)
cbar.set_ticks(np.arange(33, 37, 1))

ax4.set_ylim([100, 0])
ax4.set_xlabel('Time [h]')
ax4.set_title('Salinity [A]')

###########################################################################

# Plot 1 continued
# Second Set of Saildrones

# SST Transect SD 2
ax5 = fig.add_subplot(3, 2, 5)

cmap00 = plt.cm.jet
cmap0 = LinearSegmentedColormap.from_list('custom', cmap00(np.linspace(0.10, 1.0, 28)))
cmap = cmap0
levels = MaxNLocator(nbins=28).tick_values(10.0, 34.0)
norm = BoundaryNorm(levels, cmap0.N, clip=True)

H5 = plt.contourf(Vtime, H['Depths'], H['Verts'], levels=levels, cmap=cmap, norm=norm, extend='both')
plt.gca().invert_yaxis()

cbar = plt.colorbar(H5)
cbar.set_ticks(np.arange(12, 34, 4))

ax5.set_ylim([100, 0])
ax5.set_xlabel('Time [h]')
ax5.set_ylabel('Depth [m]')
ax5.set_title('Temperature [B]')

## SSS Transect SD 2
cmap00 = cmocean.cm.haline
cmap0 = LinearSegmentedColormap.from_list('custom', cmap00(np.linspace(0.0, 1.0, 15)))
cmap_sss = cmap0
levels = MaxNLocator(nbins=15).tick_values(32.5, 36.6)
norm_sss = BoundaryNorm(levels, cmap0.N, clip=True)

ax6 = fig.add_subplot(3, 2, 6)
H6 = plt.contourf(Vtime, I['Depths'], I['Verts'], levels=levels
                  , cmap=cmap_sss
                  , norm=norm_sss
                  , extend='both')
plt.gca().invert_yaxis()

cbar = plt.colorbar(H4)
cbar.set_ticks(np.arange(33, 37, 1))

ax6.set_ylim([100, 0])
ax6.set_xlabel('Time [h]')
ax6.set_title('Salinity [B]')

############################################################################

plt.tight_layout()
plt.savefig('Transects_Cahill', dpi=100, bbox_inches='tight')

############################################
############################################

# Plot 2 (Surface and Latent Heat Flux)

# Import data
param1 = 'swflx'
param2 = 'evap_heat'

# Read in Vertical Profiles
F = XYT_NonDepth_1VAR_trajectory(param1, ds, lon_slice, lat_slice, time_slice, trajectory)
G = XYT_NonDepth_1VAR_trajectory(param2, ds, lon_slice, lat_slice, time_slice, trajectory)
H = XYT_NonDepth_1VAR_trajectory(param1, ds, lon_slice, lat_slice, time_slice, trajectory2)
I = XYT_NonDepth_1VAR_trajectory(param2, ds, lon_slice, lat_slice, time_slice, trajectory2)

#################################################################################

# Day 1 Data

# Get data, selecting time, depth, lat/lon slice
SHflux = ds[param1].sel(time=time_slice1)
LHflux = ds[param2].sel(time=time_slice1)

# Specify longitude and latitude values for chosen domain
lons1 = SHflux.xt_ocean.values
lats1 = SHflux.yt_ocean.values

# average the data by day
avg_SHflux = SHflux.sum(dim='time') / 24  # 24 for one day
avg_LHflux = LHflux.sum(dim='time') / 24

################################################################################

# Set Up Figure
plot_area = [-115, -95, 6, 20]  # [lon1, lon2, lat1, lat2], lon is '-280'-'90' for M0M model.

cmap00 = plt.cm.jet
cmap0 = LinearSegmentedColormap.from_list('custom', cmap00(np.linspace(0.10, 1.0, 28)))
cmap = cmap0
levels = MaxNLocator(nbins=28).tick_values(-20, 320.0)
norm = BoundaryNorm(levels, cmap0.N, clip=True)

fig = plt.figure(figsize=(8.5, 7.5))

# SH Flux Map
ax1 = fig.add_subplot(2, 2, 1)
map1 = plot_map_background(plot_area=plot_area, ax=ax1)

H1 = plt.pcolormesh(lons1, lats1, avg_SHflux, cmap=cmap, norm=norm)

plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H1)
cbar.set_ticks(np.arange(0, 350, 50))

ax1.set_title('Surface Heat Flux')

# SH Flux Transect
ax3 = fig.add_subplot(2, 2, 3)

H3 = plt.plot(Vtime, F['Values'], label='Path A')
H5 = plt.plot(Vtime, H['Values'], label='Path B')

'''plt.text(0, 5, 'A/B', fontsize=18)
plt.text(34, 5, 'A\'/B\'', fontsize=18)  # Change 3/6'''

plt.legend()

# Next two lines for color limits
ax3.set_ylim([-25, 1000])
ax3.yaxis.set_ticks(np.arange(0, 1000, 200))
ax3.set_xlabel('Time [h]')
ax3.set_ylabel('Flux [$W/m^2$]')
ax3.set_title('Surface Heat Flux')

# LH Flux Map
levels2 = MaxNLocator(nbins=28).tick_values(-270, 20)
norm2 = BoundaryNorm(levels2, cmap0.N, clip=True)

ax2 = fig.add_subplot(2, 2, 2)
map2 = plot_map_background(plot_area=plot_area, ax=ax2)

H2 = plt.pcolormesh(lons1, lats1, avg_LHflux, cmap=cmap, norm=norm2)

plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H2)
cbar.set_ticks(np.arange(-250, 50, 50))

ax2.set_title('Latent Heat Flux')

# LH Flux Transect
ax4 = fig.add_subplot(2, 2, 4)
H4 = plt.plot(Vtime, G['Values'], label='Path A')
H6 = plt.plot(Vtime, I['Values'], label='Path B')

plt.legend()

ax4.set_ylim([-215, 0])
ax4.yaxis.set_ticks(np.arange(-200, 0, 50))
ax4.set_xlabel('Time [h]')
ax4.set_title('Latent Heat Flux')

##########################################################################

plt.tight_layout()
plt.savefig('Flux_Transects_Cahill', dpi=100, bbox_inches='tight')

##########################################################################
##########################################################################

# Plot 3 (Currents and Wind Speeds)

# Import data
param1 = 'u'
param2 = 'v'
param3 = 'tau_x'
param4 = 'tau_y'

# Read in Vertical Profiles
F = Depth_2VAR_trajectory(param1, param2, ds, depth, lon_slice, lat_slice, time_slice, trajectory)
G = XYU_NonDepth_2VAR_trajectory(param3, param4, ds, lon_slice, lat_slice, time_slice, trajectory)
H = Depth_2VAR_trajectory(param1, param2, ds, depth, lon_slice, lat_slice, time_slice, trajectory2)
I = XYU_NonDepth_2VAR_trajectory(param3, param4, ds, lon_slice, lat_slice, time_slice, trajectory2)


##################################################################

# Get data, selecting time, depth, lat/lon slice
datax = ds[param1].sel(time=time_slice1)
datay = ds[param2].sel(time=time_slice1)

datax = np.square(datax)
datay = np.square(datay)
data = np.add(datax, datay)
Curr = np.sqrt(data)
Curr = Curr[:, 0, :, :]

datax = ds[param3].sel(time=time_slice1)
datay = ds[param4].sel(time=time_slice1)

datax = np.square(datax)
datay = np.square(datay)
data = np.add(datax, datay)
WS = np.sqrt(data)

# Specify longitude and latitude values for chosen domain
lons1 = Curr.xu_ocean.values
lats1 = Curr.yu_ocean.values

# average the data by day
avg_Curr = Curr.sum(dim='time') / 24  # 24 for one day
avg_WS = WS.sum(dim='time') / 24

################################################################################

# Set Up Figure
plot_area = [-115, -95, 6, 20]  # [lon1, lon2, lat1, lat2], lon is '-280'-'90' for M0M model.

cmap00 = plt.cm.jet
cmap0 = LinearSegmentedColormap.from_list('custom', cmap00(np.linspace(0.10, 1.0, 28)))
cmap = cmap0
levels = MaxNLocator(nbins=28).tick_values(-.1, .8)
norm = BoundaryNorm(levels, cmap0.N, clip=True)

fig = plt.figure(figsize=(8.5, 7.5))

# Current Map
ax1 = fig.add_subplot(2, 2, 1)
map1 = plot_map_background(plot_area=plot_area, ax=ax1)

H1 = plt.pcolormesh(lons1, lats1, avg_Curr, cmap=cmap, norm=norm)

plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H1)
cbar.set_ticks(np.arange(0, .8, .1))

ax1.set_title('Surface Current Speed')

# Current Transect
ax3 = fig.add_subplot(2, 2, 3)

H3 = plt.plot(Vtime, F['Values'], label='Path A')
H5 = plt.plot(Vtime, H['Values'], label='Path B')

plt.legend()

ax3.set_ylim([-.05, .8])
ax3.yaxis.set_ticks(np.arange(0, .8, .1))
ax3.set_xlabel('Time [h]')
ax3.set_ylabel('Speed [$m/s$]')
ax3.set_title('Surface Current Speed')

# WS Map
levels2 = MaxNLocator(nbins=28).tick_values(-.05, .25)
norm2 = BoundaryNorm(levels2, cmap0.N, clip=True)

ax2 = fig.add_subplot(2, 2, 2)
map2 = plot_map_background(plot_area=plot_area, ax=ax2)

H2 = plt.pcolormesh(lons1, lats1, avg_WS, cmap=cmap, norm=norm2)

plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H2)
cbar.set_ticks(np.arange(0, .25, .05))

ax2.set_title('Wind Stress')

## WS Transect
ax4 = fig.add_subplot(2, 2, 4)
H4 = plt.plot(Vtime, G['Values'], label='Path A')
H6 = plt.plot(Vtime, I['Values'], label='Path B')

plt.legend()

ax4.set_ylim([-.01, .14])
ax4.yaxis.set_ticks(np.arange(0, .13, .02))
ax4.set_ylabel('Stress [$N/m^2$]')
ax4.set_xlabel('Time [h]')
ax4.set_title('Wind Stress')

##########################################################################

plt.tight_layout()
plt.savefig('Speeds_Transects_Cahill', dpi=100, bbox_inches='tight')

##########################################################################
##########################################################################

# Plot 4: Precip and Ocean Heights

# Import data
param1 = 'lprec'
param2 = 'eta_t'

# Read in Vertical Profiles
F = XYT_NonDepth_1VAR_trajectory(param1, ds, lon_slice, lat_slice, time_slice, trajectory)
G = XYT_NonDepth_1VAR_trajectory(param2, ds, lon_slice, lat_slice, time_slice, trajectory)
H = XYT_NonDepth_1VAR_trajectory(param1, ds, lon_slice, lat_slice, time_slice, trajectory2)
I = XYT_NonDepth_1VAR_trajectory(param2, ds, lon_slice, lat_slice, time_slice, trajectory2)

# Get data, selecting time, depth, lat/lon slice
precip = ds[param1].sel(time=time_slice1)
height = ds[param2].sel(time=time_slice1)

# Specify longitude and latitude values for chosen domain
lons1 = precip.xt_ocean.values
lats1 = precip.yt_ocean.values

# average the data by day
avg_Curr = precip.sum(dim='time') / 24  # 24 for one day
avg_WS = height.sum(dim='time') / 24

################################################################################

# Set Up Figure
plot_area = [-115, -95, 6, 20]  # [lon1, lon2, lat1, lat2], lon is '-280'-'90' for M0M model.

cmap00 = plt.cm.jet
cmap0 = LinearSegmentedColormap.from_list('custom', cmap00(np.linspace(0.10, 1.0, 28)))
cmap = cmap0
levels = MaxNLocator(nbins=28).tick_values(-.0005, .0025)
norm = BoundaryNorm(levels, cmap0.N, clip=True)

fig = plt.figure(figsize=(8.5, 7.5))

# Precip Map
ax1 = fig.add_subplot(2, 2, 1)
map1 = plot_map_background(plot_area=plot_area, ax=ax1)

H1 = plt.pcolormesh(lons1, lats1, avg_Curr, cmap=cmap, norm=norm)

plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H1)
cbar.set_ticks(np.arange(0, .0025, .0005))

ax1.set_title('Daily Precipitation')

# Precip Transect
ax3 = fig.add_subplot(2, 2, 3)

H3 = plt.plot(Vtime, F['Values'], label='Path A')
H5 = plt.plot(Vtime, H['Values'], label='Path B')

plt.legend()

ax3.set_ylim([-.00005, .0008])
ax3.yaxis.set_ticks(np.arange(0, .0008, .0001))
ax3.set_xlabel('Time [h]')
ax3.set_ylabel('Rainfall [$kg/m^2s$]')
ax3.set_title('Precipitation')

# Height Map
levels2 = MaxNLocator(nbins=28).tick_values(-.1, .6)
norm2 = BoundaryNorm(levels2, cmap0.N, clip=True)

ax2 = fig.add_subplot(2, 2, 2)
map2 = plot_map_background(plot_area=plot_area, ax=ax2)

H2 = plt.pcolormesh(lons1, lats1, avg_WS, cmap=cmap, norm=norm2)

plt.plot(trajectory['lon'], trajectory['lat'], color='gray', linewidth=2.0)
plt.text(trajectory['lon'][0], trajectory['lat'][0], 'A', color='gray')
plt.text(trajectory['lon'][143], trajectory['lat'][143], 'A\'', color='gray')

plt.plot(trajectory2['lon'], trajectory2['lat'], 'k', linewidth=2.0)
plt.text(trajectory2['lon'][0], trajectory2['lat'][0], 'B', color='k')
plt.text(trajectory2['lon'][143], trajectory2['lat'][143], 'B\'', color='k')

cbar = plt.colorbar(H2)
cbar.set_ticks(np.arange(0, .6, .1))

ax2.set_title('Ocean Surface Height')

# Height Transect
ax4 = fig.add_subplot(2, 2, 4)
H4 = plt.plot(Vtime, G['Values'], label='Path A')
H6 = plt.plot(Vtime, I['Values'], label='Path B')

plt.legend()

ax4.set_ylim([.065, .3])
ax4.yaxis.set_ticks(np.arange(.075, .3, .025))
ax4.set_ylabel('Height [m]')
ax4.set_xlabel('Time [h]')
ax4.set_title('Ocean Surface Height')

##########################################################################

plt.tight_layout()
plt.savefig('Water_Transects_Cahill', dpi=100, bbox_inches='tight')
