using Printf
using NCDatasets
using DataStructures
using CFTime, Dates
import  Statistics

function str2var(str::String, var::Any)
    str = Symbol(str)
    @eval(($str) = ($var))
end


# Define the get_gcm_info function
"""
    get_gcm_info(...)

For a specific global site, establish and store the GCM state
for each available vertical level. refers to the integer
index of the specific global site that we are interested in.
"""
# const forcingfile = "HadGEM2-A_amip.2004-2008.07"
function get_gcm_info(;
                      data_path = "/home/jbenjami/Research_Schneider/Data/cfsites/CMIP5/CFMIP2/forcing/",  
                      model     = "HadGEM2-A",
                      exper     = "amip",
                      rip       = "r1i1p1",
                      sites     = "all",
                      years     = "all",
                      months    = "all"
                      )

    @printf("--------------------------------------------------\n")
    @info @sprintf("""\n
     Experiment: GCM-driven LES(ClimateMachine)
     """)

    @printf("\n")
    print("LES = %s\n" * model * "-"* exper * "-" *rip * "-site: " * string(sites) * " years: " * string(years) * " months: " * string(months))
    @printf("--------------------------------------------------\n")
    filepath = data_path * "/" * model * "/" * exper * "/"
        # "/central/groups/esm/zhaoyi/GCMForcedLES/forcing/clima/" *
        # forcingfile *
        # ".nc"
    print(filepath)
    filenames = filter(contains(r".nc"), readdir(filepath,join=true))

    req_varnames = (
        "zg",
        "ta",
        "hus",
        "ua",
        "va",
        "pfull",
        "tntha",
        "tntva",
        "tntr",
        "tnhusha",
        "tnhusva",
        "wap",
        "hfls",
        "hfss",
	    "ts",
        "alpha",
    )
    # Load NETCDF dataset (HadGEM2-A information)
    # Load the NCDataset (currently we assume all time-stamps are 
    # in the same NCData file). We store this information in `data`. 

    # NCDataset below chokes on merging empty w/ nonempty files so we try to fix that
    # data     = Array{Any}(undef,length(filenames))
    times    = Array{Any}(undef,length(filenames)) # loads times incorrectly (seems to copy from first file)
    indices  = collect(1:length(filenames))

    for (index, file) in enumerate(filenames)
        data = NCDataset(file)
        if length(data["time"]) == 0
            indices = filter!(e->e≠index,indices) # drop those empty files
        else
            times[index] = data["time"] # save the right times
        end
    end
    
    times = times[indices]
    filenames = filenames[indices]
    time = cat(times...;dims=1) # unpacks using ellipsis

    data = NCDataset(filenames, aggdim = "time", deferopen = false)  # loads times incorrectly (seems to copy from first file)

    if months === "all" # double or triple equals? also is there an `is` operator?
        months =  Array(1:1:12)
    end
    #
    if years === "all"
        years = Dates.year.(time) # get all the years
        years = Array(minimum(years):1:maximum(years)) # get an individula vector
    end  
    time_mask = (Dates.month.(time) .∈  Ref(months)) .&  (Dates.year.(time) .∈  Ref(years)) 

    if sites === "all"
        sites = data["site"][:]
    end
    site_mask = data["site"] .∈ Ref(sites)


    # print(data)

    # To assist the user / inform them of the data processing step
    # we print out some useful information, such as groupnames 
    # and a list of available variables
    @printf("Storing information for group %s ...\n", sites)
    for (varname, var) in data #.group[groupid]
        # print(var)
        var = var=var[:,:,:] # Loading in advance makes it much faster
        for reqvar in req_varnames
            if reqvar == varname
                print("handling " * varname,"\n")
                # print(size(var),"\n")
                var = var[:,site_mask,time_mask] # seems slow for some reason, also check order
                # print(size(var),"\n")
                # Get average over time dimension
                var = Statistics.mean(var, dims = [2, 3])
                if varname == "hfls" || varname == "hfss" || varname == "ts" # surface properties, have no lev as dim 1
                    var = Statistics.mean(var, dims = 2)[2] # ???
                end
                # Assign the variable values to the appropriate converted string
                str2var(varname, var)
            end
        end
        # Store key variables
    end
    @printf("Complete\n")
    @printf("--------------------------------------------------\n")
    # @printf("Group data storage complete\n")
    return (
        zg,
        ta,
        hus,
        ua,
        va,
        pfull,
        tntha,
        tntva,
        tntr,
        tnhusha,
        tnhusva,
        wap,
        hfls,
        hfss,
        ts,
        alpha,
    )

end