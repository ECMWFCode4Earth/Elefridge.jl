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
    if true
        varname = split(split(file,"_")[end],".")[1]
        varnames[i] = varname
        println("---")
        println("Reading $varname")
        X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

        if any(X .<= 0f0)
            println("Negative/zero entries found.")
            X = X[X.>0f0]  # remove negative entries
        end

        Xlin = LinQuant16Array(X)
        Hlin[i] = bitentropy(Xlin)

        Xlog = LogQuant16Array(X,:logspace)
        Hlog[i] = bitentropy(Xlog)

        println("Hlin = $(Hlin[i]) bit, Hlog = $(Hlog[i]) bit")
    end
end

save("/Users/milan/cams/entropy/gridded_linlog16_nozeroes.jld","varnames",varnames,"Hlin",Hlin,"Hlog",Hlog)
