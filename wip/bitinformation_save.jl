using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Elefridge
using JLD

path = "/Users/milan/cams/unstructured"
filelist = filter(x->endswith(x,".grib"),readdir(path))

n = length(filelist)
varnames = fill("",n)
# ndims = ["lon","lat","vert"]
ndims = ["vert","hori"]
nbits = 32

IC = fill(0.0,n,length(ndims),nbits)

for (i,file) in enumerate(filelist)
    varname = split(split(file,"_")[end],".")[1]
    varnames[i] = varname
    println("---")
    println("Reading $varname")

    # vert x lat x lon
    X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    # convert biased exponent to signed exponent
    signed_exponent!(X)

    # permute dimensions to have lon x lat x vert, lat x lon x vert, and original
    # IC[i,1,:] = bitinformation(permutedims(X,[3,2,1]))
    # IC[i,2,:] = bitinformation(permutedims(X,[2,3,1]))
    # IC[i,3,:] = bitinformation(X)

    IC[i,1,:] = bitinformation(X)
    IC[i,2,:] = bitinformation(copy(X'))

    @save "/Users/milan/cams/entropy/information_all_unstructured.jld" varnames IC
end
