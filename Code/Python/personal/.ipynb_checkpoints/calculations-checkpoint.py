"""


"""


## Write a groupby fcn that uses bins

import numpy as np
import pandas as pd
import scipy.spatial
import itertools

import personal.math


def _structure_function_diff_power(x,inds,order): # this can work on a single pair or multidimensionally
    # should, if inds is pairs [points, pairs, other data dims] -> [points, 1, other data dims]
    
    out    = np.abs(np.diff(x.take(inds,axis=0) , axis=1))**order
    trailing_axes=tuple(np.arange(2,out.ndim,dtype=int)) # +1 so the last possible value is x.ndim which fits with our having added a dim
    # should then lead to [points, 1] bc we don't have keepdims and so lose trailing dims (note numpy needs ver 1.19.0 for keepdims)
    weight = personal.data_structures.count_value(out,axis=trailing_axes,id_fcn = lambda x: ~np.isnan(x))
    out    = np.nanmean(out, axis=trailing_axes)
    return  np.hstack((out,weight)) # shape is now hstack -> [points,2]


def structure_function(data,
                       coords,
                       random_subset_size = None, # select a random subset of the data... for speed because pairing is factorial
                       dim=None,
                       order=2,
                       groups         = None, # specify if using known_mapping or groupby_fcn, done this way to provide clarity on empty
                       mapping        = None,
                       groupby_fcn    = None, # maps data at distances to groups
                       distance_fcn   = lambda x: scipy.spatial.distance.pdist(x,metric='euclidean')[0], # should return a float
                       value_fcn      = _structure_function_diff_power,
                       structfun_data = None,
                       sorted_output  = True,
                       vectorize      = True,
                       parallelize    = False,
                       parallel_args  = {'n_jobs':4, 'require':'sharedmem'},
                       progress_bar   = False,
                       verbose        = False
                      ):
    # parellelize only parallelizes looping over groups/iterable -- you can pass vectorized functions for calculating the distances/values for example that will speed up especially vectorization behind the scenes...
    # if not vectorized, the individual loop will not benefit much from vectorization since we'll go point by point which probably is pretty small (say a small time vector for each lat/lon datapoint/row... if not you could write a vectorized fcn again still and pass it in
    """
    Structure function of order n on scalar field u defined as SF_n(r) =<|u(x+r)−u(x)|^n>, for distance r
    ------------
    If dimension 0 of data/coords has len L, this function makes use of
    > INDEX LEXICOGRAPHICAL ORDER (i.e. [(0,0),(0,1),...,(0,n), (1,0),(1,1),...,(1,L), ..., (L,0),(L,1)...,(L,L)]) for organization
    ------------
    Inputs:
    
    data          - The input data (of some np.ndarray like type) with the coordinate dimensions for the the structure function flattened
                    over the 0'th dimension.
                    E.g. if we have data in x,y,z,t and want the structure function over x,y, we should flatten the x and y dimensions,
                    yielding an array with shape [s_0, s_1, s_2] = [s_x*s_y, s_z, s_t]
                    If necessary do this outside this function, possibly with flatten_dimensions from personal.data_structures
                    
    coords        - The input coordinates as an array of shape [s_0, n_coords] where s_0 is the                 
                    > This fcn uses matrix arithmetic to speed itself up by broadcasting arithmetic over trailing dimensions
                    
                    
    order         - The structure function order
    bins          - We can create distance groupings using bins, default (None) lists every found distance uniquely
    mapping - If we've already calculated the mappings between position in array to distance groups in our structure function output,
                    this can be used to save time.
                    You should provide simply the mapping as an array containing the mapping between groups and distances,
                    and a list of group indices to match the data pairs we have in index lexographical order 
                    e.g. with groups [  0.1,   0.33,   1.40,   3.19] mapping could be [0.1, 1.40, 0.1, 0.1, 3.19, 0.33,...]
                      or with groups ['1', '2', '3', '4']            mapping could be ['4', '2', '1', '4', '3',...]
                      
                    Note since the # of pairings could be quite large, one could imagine using a generator that yields these values
                    based on some underlying fcn
                    
                    
    groupby_fcn   - Default is None, but accepts fcn of groupby_fcn(dist)
                    Useful for example for grouping data into bins or other such collections
                    Must return group i.e. ([0.33,1.40]) or (1.40) or 2 '2' or whatever
                    
    groups        - Default will just use the set of all outputs from groupby_fcn 
                    
    distance_fcn  - Calculates the distances between all pairs of points in index lexicographic order
                    > You can pass any function as metric that calculates distances between pairs of points (i.e. one for lat/lon points)
                      I decided to keep this as a lambda fcn rather than hard coding pdist in case you want to pass kwargs to pdist
                      or use something completely different
    ------------
    Output has form {dist:value}, note dist could be a bin (x1,x2) or something like that, as long as hashable
    ------------
    
    Regardless of input size dimension, process is as follows:
    
    In INDEX LEXICOGRAPHICAL ORDER (i.e. [(0,0),(0,1),...,(0,n), (1,0),(1,1),...,(1,n), ..., (n,0),(n,1)...,(n,n)]) along dimension 0:
    > generate unique pairings itertools or something else? apply metric_fcn (i.e. abs(mean()) to get the sorting criteria
    > ... or provide 
    groupby either the bins or unique value? if you already have the ravel_indices you could have a shortcut grouping criteria
    
    For applying along a dimension, use apply_along_dimension or some similar tactic... (can speed up if the distance mapping is repeated)
    
    Remains to be seen how plays wit dask but should work? (pure dask not xarray)
    
    We assume if the full vectors don't fit in memory, and iterator/generator form is too slow to be tenable....
    If you wish this to work in this way, try memmap'd arrays, chunking, dask, precalculating the fcn, or a random point subset
    """
    

   

    
        

    # handle length, and take a random subset if necessary
    L = len(data)
    if random_subset_size is not None: # take subset
        print('taking random data subset')
        selections = np.random.choice(L, random_subset_size, replace=False)
        L    = random_subset_size
        data = data[selections]
        coords = coords[selections]
        # since mapping, structfun_data would be defined along the pairs, we will square form and drop the rows and columns not in selections, assuming index lexicographic order!!
        if mapping is not None:
            mapping        = scipy.spatial.distance.squareform(mapping)
            mapping        = mapping[selections][:,selections]
            mapping        = scipy.spatial.distance.squareform(mapping)
        if structfun_data is not None:
            structfun_data = scipy.spatial.distance.squareform(structfun_data)
            structfun_data = structfun_data[selections][:,selections]
            structfun_data = scipy.spatial.distance.squareform(structfun_data)
    
        
    # calculate pairs
    if mapping is None:
        if vectorize:
            pairs = personal.math.combinations(range(L),2,allow_repeats=False) # same as itertools.combinations(range(L),2) but full array
        else:
            pairs = itertools.combinations(range(L),2) 
#         npairs = L*(L-1)//2 # L choose 2  

    # handle known and unknown options and set up grouping and distances
    if vectorize:
        if mapping is None:
            mapping = distance_fcn(coords) # for pairs, a fcn to generate all pair distances in lexicographic order or some other mapping
            if grouby_fcn:
                mapping = groupby_fcn(mapping) # replace dists with mapping to whatever the groupby_fcn returns
            
        if structfun_data is None:
            structfun_data = value_fcn(data, personal.math.combinations(range(L),2,allow_repeats=False), order) # use on array of pairs of points
    else:
        if mapping is None:
            mapping = map (lambda x: distance_fcn(np.vstack((x[0],x[1]))),itertools.combinations(coords,2) ) # mapping to coord pair list
            if groupby_fcn:
                mapping = map(groupby_fcn,mapping)  # replace dists with mapping to whatever the groupby_fcn returns (this is a generator!)
                
        if structfun_data is None:
            structfun_data = map(lambda pair: value_fcn(data,np.array([pair]),order), itertools.combinations(range(L),2))
        
   
    if groups is None: # default to span of groupby_fcn output
        if vectorize:
            groups = np.unique(mapping)
            struct_fun_info = {key:[0,0] for key in groups}
        else:
            struct_fun_info = dict() # start empty so we can use our generators to fill
    else:
        struct_fun_info = {key:[0,0] for key in groups} # from keys would fucc up 
    
    
    if vectorize:
        it = groups # iterate over the distsance_bin groupings we calculated
        tqdm_total = len(groups)
    else:
        it = zip(mapping,structfun_data) # iterate over each result (note this uses a different parallel_fcn later) structfun_data is an array w/ rows [value,weight] where value might be nan for weight 0
        tqdm_total = L*(L-1)//2
        
    if progress_bar:
#         from tqdm import tqdm
        from tqdm.auto import tqdm
        it = tqdm(it,total=tqdm_total)

    if parallelize:      
        print('parallelizing')
        parallel_args.update({'require':'sharedmem'}) #ensure shared memory for dict
        from joblib import Parallel, delayed

    # Put it all together
    if vectorize: # we'll already know the groups...
        if parallelize:
            # write fcn to take dst_group,mapping=mapping,structfun_data=structfun_data and update struct_fun_info
            def parallel_fcn(dst_grp, struct_fun_info = struct_fun_info, structfun_data=structfun_data,mapping=mapping):
                dg =  structfun_data[mapping==dst_grp]
#                 print(dst_grp,dg)
                struct_fun_info[dst_grp] = np.nansum(np.prod(dg,axis=1)) / np.sum(dg[:,1])
                return # return nothing
            Parallel(**parallel_args)(delayed(parallel_fcn)(dst_grp,struct_fun_info = struct_fun_info, structfun_data=structfun_data,mapping=mapping) for dst_grp in it)
        else: # no parallelize
            for dst_grp in it:
                dg =  structfun_data[mapping==dst_grp]
                struct_fun_info[dst_grp] = np.nansum(np.prod(dg,axis=1)) / np.sum(dg[:,1]) # weighted mean
    else:
        if parallelize:      
            # write fcn to take dst_group,mapping=mapping,structfun_data=structfun_data and update struct_fun_info
            def parallel_fcn(dst_grp, strfcn_datum, struct_fun_info = struct_fun_info):
                new = np.vstack((strfcn_datum,struct_fun_info.get(dst_grp,[0,0])))
                new_weight =  np.sum(new[:,1])
                struct_fun_info[dst_grp] = [np.nansum(np.prod(new,axis=1))/new_weight, new_weight] # update our dict, memory mapped
                return # return nothing
                                
            Parallel(**parallel_args)(delayed(parallel_fcn)(dst_grp,strfcn_datum,struct_fun_info = struct_fun_info) for dst_grp,strfcn_datum in it)
            struct_fun_info = {key:val[0] for key,val in struct_fun_info.items()} # drop the weights at the end


        else: # no parallelize
            for dst_grp,strfcn_datum in it: # iterate over generator/iterator/iterable
                struct_fun_info[dst_grp] = weighted_nanmean(strfcn_dat,struct_fun_info.get(dst_grp,[0,0])) # update value
   
    
    
    struct_fun_info = pd.Series(struct_fun_info)
    return struct_fun_info.sort_index() if sorted_output else output



