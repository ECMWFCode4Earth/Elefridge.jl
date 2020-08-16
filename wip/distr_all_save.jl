using Statistics, StatsBase, LinearAlgebra
using PyCall
xr = pyimport("xarray")
using PyPlot
using Elefridge
using JLD

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,".grib"),readdir(path))

n = length(filelist)
varnames = fill("",n)

# example histogram for preallocation
Hex = fit(Histogram,rand(Float32,100),[0f0,0.5f0,1f0])
S = fill(normalize(Hex,mode=:density),n)

for (i,file) in enumerate(filelist)
    varname = split(split(file,"0_")[end],".")[1]
    varnames[i] = varname
    println("---")
    println("Reading $varname")
    local X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    if all(X .== zero(eltype(X)))
        println("Only zeros in $varname")
    else
        mi = log(Elefridge.minpos(X))
        ma = log(maximum(X))

        local n = 300
        Δ = (ma-mi)/n

        bins = vcat([-Inf32,0f0,floatmin(Float32)],exp.(mi:Δ:ma),[prevfloat(Inf32),Inf32])

        local H = fit(Histogram,vec(X),bins)

        # normalize by overall elements not by area
        Hn = normalize(H,mode=:density)
        Hn.weights .= H.weights / prod(size(X))
        S[i] = Hn
    end
end

@save "/Users/milan/cams/entropy/distribution_all_gridded.jld" varnames S
