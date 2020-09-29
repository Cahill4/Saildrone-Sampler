#import xarray as xr
from scipy.interpolate import NearestNDInterpolator
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
#from matplotlib.ticker import MaxNLocator

##########################################################

# Define calc_trajectory function
def calc_trajectory(speed, duration_hours, direction_ccw_from_east,
                    prev_trajectory=None, dt_start=None, lon_start=None, lat_start=None):
    ##
    ## A trajectory dictionary has keys 'datetime','lon','lat'
    ##
    factor = 3600.0 / 111000.0  # m/s --> deg / hour at equator. (For lon, adjust for latitude below.)
    direction_radians = direction_ccw_from_east * 3.14159 / 180.0

    if prev_trajectory is None:
        lat_radians = lat_start * 3.14159 / 180.0
        hours = range(int(duration_hours) + 1)
        dt_track = [dt_start + dt.timedelta(hours=x) for x in hours]
        lon_track = np.array(
            [lon_start + factor * np.cos(lat_radians) * speed * x * np.cos(direction_radians) for x in hours])
        lat_track = np.array([lat_start + factor * speed * x * np.sin(direction_radians) for x in hours])
    else:
        lat_radians = prev_trajectory['lat'][-1] * 3.14159 / 180.0
        dt_start = prev_trajectory['datetime'][-1]
        hours = range(1, int(duration_hours) + 1)
        dt_new_track = [dt_start + dt.timedelta(hours=x) for x in hours]
        dt_track = prev_trajectory['datetime'] + dt_new_track
        lon_track = np.append(prev_trajectory['lon']
                              , np.array(
                [prev_trajectory['lon'][-1] + factor * np.cos(lat_radians) * speed * x * np.cos(direction_radians) for x
                 in hours]))
        lat_track = np.append(prev_trajectory['lat']
                              , np.array(
                [prev_trajectory['lat'][-1] + factor * speed * x * np.sin(direction_radians) for x in hours]))

    ## Prepare output dictionary.
    F = {}
    F['datetime'] = dt_track
    F['lon'] = lon_track
    F['lat'] = lat_track
    return F

##########################################################

def vertical_prof_trajectory(param, ds, depth_slice, lon_slice, lat_slice, time_slice, trajectory):

    # Define data
    data = ds[param].sel(st_ocean=depth_slice, xt_ocean=lon_slice, yt_ocean=lat_slice, time=time_slice)

    ## Handle points if only one point specified.
    points = np.array([trajectory['lon'].tolist(), trajectory['lat'].tolist()]).T

    # Define lat and lons from data
    lon = data.xt_ocean.values
    lat = data.yt_ocean.values

    # Stack lat and lons for Nearest Neighbor
    stackedlat = np.repeat(lat, 40)
    stackedlon = np.tile(lon, 48)

    nz = len(data.st_ocean.values)

    field_profile_list = np.zeros([nz, points.shape[0]])

    # Process each field at each time.
    for tt in range(len(trajectory['datetime'])):
        datedata = data[tt, :, :, :]

        for kk in range(nz):
            stacked_data = datedata[kk, :, :].stack(z=("yt_ocean", "xt_ocean"))
            interp = NearestNDInterpolator((stackedlon, stackedlat), stacked_data)
            field_profile_list[kk, tt] = interp(points[tt, :])

    depths = data.st_ocean.values

    # Prepare Outputs
    FOUT = {}
    FOUT['Verts'] = field_profile_list
    FOUT['Depths'] = depths

    return FOUT

# Xu/Yu (extra lat/lon)
def Depth_2VAR_trajectory(param1, param2, ds, depth, lon_slice, lat_slice, time_slice, trajectory):

    # Define data
    data1 = ds[param1].sel(xu_ocean=lon_slice, yu_ocean=lat_slice, time=time_slice)
    data2 = ds[param2].sel(xu_ocean=lon_slice, yu_ocean=lat_slice, time=time_slice)

    # Find magnitude of Windstress
    data1 = np.square(data1)
    data2 = np.square(data2)
    data = np.add(data1, data2)
    data = np.sqrt(data)

    ## Handle points if only one point specified.
    points = np.array([trajectory['lon'].tolist(), trajectory['lat'].tolist()]).T

    # Define lat and lons from data
    lon = data.xu_ocean.values
    lat = data.yu_ocean.values

    # Stack lat and lons for Nearest Neighbor
    stackedlat = np.repeat(lat, 41)
    stackedlon = np.tile(lon, 49)

    field_profile_list = np.zeros([points.shape[0]])

    # Process each field at each time.
    for tt in range(len(trajectory['datetime'])):
        stacked_data = data[tt, depth, :, :].stack(z=("yu_ocean", "xu_ocean"))
        interp = NearestNDInterpolator((stackedlon, stackedlat), stacked_data)
        field_profile_list[tt] = interp(points[tt, :])

       # Prepare Outputs
    FOUT = {}
    FOUT['Values'] = field_profile_list

    return FOUT

def XYU_NonDepth_2VAR_trajectory(param1, param2, ds, lon_slice, lat_slice, time_slice, trajectory):

    # Define data
    data1 = ds[param1].sel(xu_ocean=lon_slice, yu_ocean=lat_slice, time=time_slice)
    data2 = ds[param2].sel(xu_ocean=lon_slice, yu_ocean=lat_slice, time=time_slice)

    # Find magnitude of Windstress
    data1 = np.square(data1)
    data2 = np.square(data2)
    data = np.add(data1, data2)
    data = np.sqrt(data)

    ## Handle points if only one point specified.
    points = np.array([trajectory['lon'].tolist(), trajectory['lat'].tolist()]).T

    # Define lat and lons from data
    lon = data.xu_ocean.values
    lat = data.yu_ocean.values

    # Stack lat and lons for Nearest Neighbor
    stackedlat = np.repeat(lat, 41)
    stackedlon = np.tile(lon, 49)

    field_profile_list = np.zeros([points.shape[0]])

    # Process each field at each time.
    for tt in range(len(trajectory['datetime'])):
        stacked_data = data[tt, :, :].stack(z=("yu_ocean", "xu_ocean"))
        interp = NearestNDInterpolator((stackedlon, stackedlat), stacked_data)
        field_profile_list[tt] = interp(points[tt, :])


    # Prepare Outputs
    FOUT = {}
    FOUT['Values'] = field_profile_list

    return FOUT

def XYT_NonDepth_1VAR_trajectory(param, ds, lon_slice, lat_slice, time_slice, trajectory):

    # Define data
    data = ds[param].sel(xt_ocean=lon_slice, yt_ocean=lat_slice, time=time_slice)

    ## Handle points if only one point specified.
    points = np.array([trajectory['lon'].tolist(), trajectory['lat'].tolist()]).T

    # Define lat and lons from data
    lon = data.xt_ocean.values
    lat = data.yt_ocean.values

    # Stack lat and lons for Nearest Neighbor
    stackedlat = np.repeat(lat, 40)
    stackedlon = np.tile(lon, 48)

    field_profile_list = np.zeros([points.shape[0]])

    # Process each field at each time.
    for tt in range(len(trajectory['datetime'])):
        stacked_data = data[tt, :, :].stack(z=("yt_ocean", "xt_ocean"))
        interp = NearestNDInterpolator((stackedlon, stackedlat), stacked_data)
        field_profile_list[tt] = interp(points[tt, :])

    # Prepare Outputs
    FOUT = {}
    FOUT['Values'] = field_profile_list

    return FOUT

def plot_map_background(plot_area=[-180, 180, -70, 70], ax=plt.gca()
                        , projection='cyl', resolution='l'
                        , draw_meridians=np.arange(0, 360, 10), meridians_line_width=0.25
                        , meridians_labels=[0, 0, 0, 1], meridians_fontsize=10
                        , meridians_dashes=[4, 2]
                        , draw_parallels=np.arange(-90, 90, 10), parallels_line_width=0.25
                        , parallels_labels=[1, 0, 0, 0], parallels_fontsize=10
                        , parallels_dashes=[4, 2]
                        , coastline_line_width=0.5
                        , countries_line_width=0.25
                        , x_tick_rotation=45.0
                        , y_tick_rotation=45.0):
    map1 = Basemap(projection=projection,
                   llcrnrlon=plot_area[0],
                   urcrnrlon=plot_area[1],
                   llcrnrlat=plot_area[2],
                   urcrnrlat=plot_area[3],
                   resolution=resolution, ax=ax)

    # draw lat/lon grid lines degrees.
    meridians = map1.drawmeridians(draw_meridians, linewidth=meridians_line_width
                                   , labels=meridians_labels, fontsize=meridians_fontsize
                                   , dashes=meridians_dashes, textcolor='white')
    parallels = map1.drawparallels(draw_parallels, linewidth=parallels_line_width
                                   , labels=parallels_labels, fontsize=parallels_fontsize
                                   , dashes=parallels_dashes, textcolor='white')
    map1.drawcoastlines(linewidth=coastline_line_width)
    map1.drawcountries(linewidth=countries_line_width)
    map1.fillcontinents(color='white', lake_color='black')

    for m in meridians:
        try:
            meridians[m][1][0].set_rotation(x_tick_rotation)
        except:
            pass

    for p in parallels:
        try:
            parallels[p][1][0].set_rotation(y_tick_rotation)
        except:
            pass

    return map1
