import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from shapely.geometry.polygon import Polygon
import numpy as np
import os
import xarray as xr
import dask
import dask.array as da
import scipy
import sklearn

import personal.data_structures
import personal.math
import personal.calculations
import personal.system
import personal.plots
import personal.constants

#%% Plot setup
dpi = 100
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

coords = {'domain_lons': [ x-360 for x in [234.3, 237.1, 237.6, 234.9]], # [-125.7, -122.9, -122.4, -125.1]
          'domain_lats': [035.5, 038.0, 037.2, 034.7]
         }

def draw_SMODE_box(ax):
    """
    Draws the S-MODE box over an existing area (maybe overload with a general draw box fcn?
    """
    # top left, top right, bottom right, bottom left
    domain_lons = [234.3, 237.1, 237.6, 234.9]
    domain_lons = [x-360 for x in domain_lons]
    domain_lats = [035.5, 038.0, 037.2, 034.7]
    domain_lats.append(domain_lats[0])
    smode_box = (*zip(domain_lons, domain_lats),)
    pgon = Polygon(smode_box)
    return ax.add_geometries([pgon], crs=ccrs.PlateCarree(), facecolor='None', edgecolor='black', linestyle='dashed', linewidth=1.2,alpha=1.0)


def select_SMODE_region(data,rechunk=True,time_chunk=1):
    """
    Selects data in a subset box region that contains S-MODE area (maybe overload with a general select box fcn?)
    """
    
    lats = [34,39]
    lons = [235.5, 239.5]; lons = [x-360 for x in lons]
    # Return normal data slices and rechunk to be 1 in time and lat/lon slices
    out = data.sel(lat=slice(*lats),lon=slice(*lons))
    return out.chunk({'time':time_chunk,'lat':len(out['lat']),'lon':len(out['lon'])}) if rechunk else out

def set_ax_SMODE_region(ax):
    ax.set_ylim([36, 39]);
    ax.set_xlim(np.array([235.5,239.5])-360)
    return
    
def SMODE_xr_plot(xr_data,dpi=dpi,**kwargs):
    """
    Plots in the SMODE region provided data is xr_data form lat x lon only
    """
    import matplotlib.colors
    from matplotlib import cm
    
    fig = plt.figure(figsize=(960/dpi,900/dpi),dpi=dpi)

    ax = plt.axes(projection=ccrs.PlateCarree(),facecolor='black')
    ax.coastlines(resolution='10m', color='black', linewidth=1);
    land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='b',facecolor='black')
    ax.background_patch.set_fill(False)
    ax.add_feature(land_10m, zorder=100, edgecolor='k')
    g =  ax.gridlines(draw_labels=True)
    g.xlabels_top=False
    g.ylabels_right=False

    draw_SMODE_box(ax)
    set_ax_SMODE_region(ax)
    plt.rcParams.update({'font.size': 12})
    plt.subplots_adjust(bottom=0.05, top=.95, left=0.05, right=0.99)

    cax = xr_data.plot(ax=ax,**kwargs)

    return fig,cax

# To Do -- add version for just regular data



def calculate_structure_fcn(
    data                    = None,              # dictionary of form {'data':data,'precomputed':bool}
    precomputed             = {},   # Things to include include coords, distances, bin_mapping
    precomputed_mask        = None, # can be a string pointing out a precomputed mask file...
    save_filepath           = None,
    save_variables          = {'data':None, 'mask':None,'distances':None,'bin_mapping':None, 'structure_fcn':None}, # save items in this dictionary to given filenames
    time_groups             = {'all':slice(0,None,None)}, # how to group the time, can be a list of logical arrays
    metadata                = None,
    nbins                   = 50,
    order                   = 2,
    random_subset_size      = None,
    fill_output_na          = False,
    file_path               = os.path.abspath(os.path.dirname('')),
    data_relpath            = "../Data/HF_Radar/2km/processed/smode_region/",
    output_savepath         = "../Data/HF_Radar/2km/processed/smode_region/miscellaneous",
    ch_size_fcn             = None,
    vectorize               = True,
    parallelize             = True,
    parallel_args           = {'n_jobs':20, 'require':'sharedmem'},
    progress_bar            = True,
    value_fcn_parallelize   = True,
    value_fcn_parallel_args = None,
    value_fcn_progress_bar  = True
    ):
    
    print('precomputing...')
    out = process_for_structure_fcn(data               = data,
                                    precomputed        = precomputed,
                                    save_filepath      = None, # don't save her
                                    save_variables     = None, # don't save here
                                    metadata           = metadata,
                                    nbins              = nbins,
                                    random_subset_size = random_subset_size,
                                    file_path          = file_path, 
                                    data_relpath       = data_relpath,
                                    output_savepath    = output_savepath
                                    ) # will return a dict with everything in it....
    random_subset_size = None # should have been already processed
    
    a = personal.constants.earth['radius']    
    data                = out['data']  
    mask                = out['mask']
    distances           = out['distances']  
    distance_bins       = out['distance_bins']  
    distance_bin_ranges = out['distance_bin_ranges']  
    bin_mapping         = out['bin_mapping']  
    coords              = out['coords']  
    
    

        
        


    print('calculating structure function...')        
    # Helper fcn we want so we want, for structure_fcn
    if ch_size_fcn is None: ch_size_fcn = lambda data: round(6000 * (max(1,70000/data.shape[1]))**1.0)
        
        
    def value_fcn(data,inds,order,ch_size_fcn = ch_size_fcn,parallelize=value_fcn_parallelize,parallel_args=value_fcn_parallel_args,progress_bar=value_fcn_progress_bar ): # ch_size_fcn tuned based on experience for eady...
        def diff_power(inds,out,order): # this can work on a single pair or multidimensionally, argument order changed so inds can be first for parallelization
            # should, if inds is pairs [points, pairs, other data dims] -> [points, 1, other data dims]
#             out    = np.abs(np.diff(out.take(inds,axis=0) , axis=1))**order 
            out    = np.abs(np.diff(out[inds] , axis=1))**order # Use fancy indexing, then diff along new axis, then abs and raies to power
            trailing_axes=tuple(np.arange(2,out.ndim,dtype=int))
            # should then lead to [points, 1] bc we don't have keepdims and so lose trailing dims (note numpy needs ver 1.19.0 for keepdims)
            weight = personal.data_structures.count_value(out,axis=trailing_axes,id_fcn = lambda x: ~np.isnan(x)) # collapses trailing axes...
            out    = np.nanmean(out, axis=trailing_axes)
            return  np.hstack((out,weight)) # shape is now hstack -> [points,2]
        ch_size = ch_size_fcn(data)
        if parallel_args is None: parallel_args={'n_jobs':min(20,int(np.ceil(len(inds)/ch_size)))}
        print('chunk size: ' + str(ch_size) + ' for total length ' + str(len(inds)) + ' -> ' + str(np.ceil(len(inds)/ch_size)) + ' chunks')
        print(parallel_args)
        return personal.data_structures.apply_in_chunks(diff_power, inds, data, order,
                                                        output_shape=(len(inds),2), input_chunk_sizes = {0:ch_size}, output_chunk_sizes={},
                                                        parallelize=parallelize, parallel_args=parallel_args, progress_bar=progress_bar ) # fill in remaining args wit *args (see apply_in_chunks documentation)

    structure_fcn = dict.fromkeys(time_groups.keys())
    for time_group, indices in time_groups.items():
        print('Processing: ' + str(time_group))
        n = personal.calculations.structure_function(data[:,indices],coords,
                                                     mapping = bin_mapping, groups = list(range(nbins)),random_subset_size=random_subset_size,
                                                     parallelize=parallelize,vectorize=vectorize,progress_bar=progress_bar,value_fcn=value_fcn, parallel_args  = parallel_args)
        
        n.index = personal.math.geometric_mean(distance_bin_ranges,axis=1) # geometric mean of bin edges
        if fill_output_na:
            structure_fcn[time_group] = n.fillna(method='ffill')
        else:
            structure_fcn[time_group] = n
    
    out.update({'structure_fcn':structure_fcn})
#     out = {'structure_fcn':struture_fcn, 'data':data, 'mask':mask,'distances':distances,'distance_bins':distance_bins, 'distance_bin_ranges':distance_bin_ranges,'bin_mapping':bin_mapping,'coords':coords}
    if save_filepath is not None:
        print('saving...')        
        personal.IO.pickle_save( personal.data_structures.subset_dict(out,save_variables) ,save_filepath)
        


    return out


def process_for_structure_fcn(
    data               = None,              # the data to pre-process (should be an xarray)
    precomputed        = {},   # Things to include include coords, distances, bin_mapping (also if you dont want something computed, add it to precomputed, but remove from save_variables)
    precomputed_mask   = None,    # can be a string pointing out a precomputed mask file...
    save_filepath      = None,
    save_variables     = {'data':None, 'mask':None,'distances':None,'bin_mapping':None,'coords':None}, # save items in this dictionary to given filenames
    time_groups        = {'all':slice(0,None,None)}, # how to group the time, can be a list of logical arrays
    metadata           = None,
    nbins              = 50,
    random_subset_size = None,
    file_path          = os.path.abspath(os.path.dirname('')), # not sure why it's returning tuple inside function, take first index
    data_relpath       = "../Data/HF_Radar/2km/processed/smode_region/",
    output_savepath    = "../Data/HF_Radar/2km/processed/smode_region/miscellaneous"
    ):
    
    
    # Handle d (the only thing that doesn't have to be precomputed...)
    mask     = precomputed.get('mask',None)
    if precomputed.get('data',None) is None:
        d_type = type(data)
        if d_type is xr.core.dataarray.DataArray: # is xarray DataArray
            metadata = data # use metadata from a fully operational DataArray
            print('getting values from DataArray')
            data = data.values # .data would return dask or numpy array
            print('Done')
        elif   d_type is da.core.Array: # is dask array
            data = data.compute()
        elif d_type is np.ndarray: # is numpy array
            pass
        else:
            raise TypeError('Unsupported data type ' + str(d_type))
        
        print('transposing')
        data = data.transpose([1,2,0]) # make dims be lat,lon,time
        print('reshaping')
        data = data.reshape(-1, *data.shape[-1:]) # then flatten to list along coords (flatten all but last dim time)
        if mask is None: mask =  ~np.all(np.isnan(data), axis=1)
        test_mask = ~np.all(data==0, axis=1)
        data = data[mask]
    else:
        data = precomputed['data']
            
      
    a = personal.constants.earth['radius']    
    
    
    if metadata is None:
        metadata = xr.open_mfdataset(os.path.normpath(os.path.join(file_path, output_savepath)) + '/HFRADAR_SMODE_Region_2km_Resolution_Hourly_Metadata.nc', decode_cf=True, decode_times=True,combine='by_coords',data_vars="minimal", coords="minimal", compat="override") # chunk doesnt expand beyond file boundaries
            
    if precomputed.get('coords', None) is None:
        print('calculating coordinates')
        lat  = metadata['lat']
        lon  = metadata['lon']
        latg,long = np.meshgrid(lat,lon,indexing='ij') # 'ij' to match the shape of x,y
        coords = np.float32(np.stack((latg.flatten(), long.flatten()),axis=1)) # check order (x,y matches lon,lat)
        coords = coords[mask] 
    else:
        coords = precomputed['coords']
        
        
    L = len(data)
    if random_subset_size is not None: # take subset (asusumes anything processed has all the nan masking applied etc)
        print('taking data subset')
        selections = np.random.choice(L, random_subset_size, replace=False)
        L          = random_subset_size # replace old value
        data       = data[selections]
        coords     = coords[selections]
        # since mapping, structfun_data would be defined along the pairs, we will square form and drop the rows and columns not in selections, assuming index lexicographic order!!
        if precomputed.get('distances', None) is not None:
            pass # handled later          
        if precomputed.get('mapping', None) is not None:
            pass # handled later
        if precomputed.get('structfun_data', None) is not None:
            structfun_data = scipy.spatial.distance.squareform(structfun_data)
            structfun_data = structfun_data[selections][:,selections]
            structfun_data = scipy.spatial.distance.squareform(structfun_data)

        
    # Distances and Bins
    if precomputed.get('distances',None) is None:
        print('calculating distances')
        distances = scipy.spatial.distance.squareform(sklearn.metrics.pairwise.haversine_distances(np.radians(coords)), checks=False)*a # calc can run out of memory, but isn't as slow... consider a cython implementation, flip for haversine
    else:
        distances = precomputed['distances']
        if random_subset_size is not None: # take subset (asusumes anything processed has all the nan masking applied etc)
            distances      = scipy.spatial.distance.squareform(distances)
            distances      = distances[selections][:,selections]
            distances      = scipy.spatial.distance.squareform(distances)  

        
    if precomputed.get('distance_bins',None) is None:
        print('binning')
        # furthest_distance = personal.geometry.haversine_np([long[0][0],latg[0][0]],[long[-1][-1],latg[-1][-1]],a=a);
        furthest_distance   = np.max(distances)
        shortest_distance   = np.min(distances) 
        distance_bins       = np.logspace(np.log10(shortest_distance-1e-10), np.log10(furthest_distance),num=nbins+1,endpoint=True)
    else:
        distance_bins       = precomputed['distance_bins']
        
    if precomputed.get('distance_bin_ranges',None) is None:
        distance_bin_ranges = np.vstack((distance_bins[:-1],distance_bins[1:])).T 
    else:
        distance_bin_ranges = precomputed['distance_bin_ranges']
        
    if precomputed.get('bin_mapping',None) is None:
        bin_mapping = np.digitize(distances,distance_bins,right=True) - 1 #  -1 bc digitize only goes down to 1 and we want bins to match up with indices of bin_ranges
    else:
        bin_mapping = precomputed['bin_mapping']
        if random_subset_size is not None: # take subset (asusumes anything processed has all the nan masking applied etc)
            bin_mapping        = scipy.spatial.distance.squareform(bin_mapping)
            bin_mapping        = bin_mapping[selections][:,selections]
            bin_mapping        = scipy.spatial.distance.squareform(bin_mapping)
        
            
    out = {'data':data, 'mask':mask,'distances':distances,'distance_bins':distance_bins, 'distance_bin_ranges':distance_bin_ranges,'bin_mapping':bin_mapping,'coords':coords}
    
    if save_filepath is not None:
        personal.IO.pickle_save( personal.data_structures.subset_dict(out,save_variables) ,save_filepath)
        
    return out

    
    
