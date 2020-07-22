using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Elefridge
using JLD
using BFloat16s

path = "/Users/milan/cams/gridded"
filelist = filter(x->endswith(x,".grib"),readdir(path))

n = length(filelist)
varnames = Array{String,1}(undef,n)
Hlin = Array{Float64,1}(undef,n)
Hlog = Array{Float64,1}(undef,n)

for (i,file) in enumerate(filelist)
    varname = split(split(file,"_")[end],".")[1]
    varnames[i] = varname
    println("---")
    println("Reading $varname")
    X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    try
        Xlin = LinQuant16Array(X)
        Hlin[i] = bitentropy(Xlin)
    catch
        Hlin[i] = 0
    end

    try
        Xlog = LogQuant16Array(X)
        Hlog[i] = bitentropy(Xlog)
    catch
        Hlog[i] = 0
    end

    println("Hlin = $(Hlin[i]) bit, Hlog = $(Hlog[i]) bit")
end

save("/Users/milan/cams/entropy/gridded_linlog16.jld","varnames",varnames,"Hlin",Hlin,"Hlog",Hlog)
