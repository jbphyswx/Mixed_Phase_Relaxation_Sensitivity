# We first specify the NetCDF file from which we wish to read our 
# GCM values.
# Utility function to read and store variables directly from the 
# NetCDF file
"""
    str2var(str::String, var::Any)
Helper function allowing variables read in from the GCM file
to be made available to the LES simulation.
"""
function str2var(str::String, var::Any)
    str = Symbol(str)
    @eval(($str) = ($var))
end

# Define the get_gcm_info function
"""
    get_gcm_info(groupid)
For a specific global site, establish and store the GCM state
for each available vertical level. `groupid` refers to the integer
index of the specific global site that we are interested in.
"""
const forcingfile = "HadGEM2-A_amip.2004-2008.07"
function get_gcm_info(groupid)

    @printf("--------------------------------------------------\n")
    @info @sprintf("""\n
     Experiment: GCM(HadGEM2-A) driven LES(ClimateMachine)
     """)

    @printf("\n")
    @printf("HadGEM2-A_LES = %s\n", groupid)
    @printf("--------------------------------------------------\n")
    filename =
        "/central/groups/esm/zhaoyi/GCMForcedLES/forcing/clima/" *
        forcingfile *
        ".nc"

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
    data = NCDataset(filename)
    # To assist the user / inform them of the data processing step
    # we print out some useful information, such as groupnames 
    # and a list of available variables
    @printf("Storing information for group %s ...", groupid)
    for (varname, var) in data.group[groupid]
        for reqvar in req_varnames
            if reqvar == varname
                # Get average over time dimension
                var = mean(var, dims = 2)
                if varname == "hfls" || varname == "hfss" || varname == "ts"
                    var = mean(var, dims = 1)[1]
                end
                # Assign the variable values to the appropriate converted string
                str2var(varname, var)
            end
        end
        # Store key variables
    end
    @printf("Complete\n")
    @printf("--------------------------------------------------\n")
    @printf("Group data storage complete\n")
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
