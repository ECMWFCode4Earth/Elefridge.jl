using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using PyPlot
using Elefridge
using JLD
using BFloat16s

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,".grib"),readdir(path))

n = length(filelist)
varnames = fill("",n)

for (i,file) in enumerate(filelist)
    varname = split(split(file,"0_")[end],".")[1]
    varnames[i] = varname
    println("---")
    println("Reading $varname")
    X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    N = prod(size(X))
    nzero = sum(X .== zero(eltype(X)))
    nneg = sum(X .< zero(eltype(X)))

    # println("$nzero, $(nzero/N*100)% 0 entries.")
    # println("$nneg, $(nneg/N*100)% negative entries.")

    H = bitentropy(X)
    println("$H bit Float32 entropy.")
end
