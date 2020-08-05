using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Blosc
using Elefridge
using JLD
using ZfpCompression

path = "/Users/milan/cams/gridded"
filelist = filter(x->endswith(x,".grib"),readdir(path))

n = length(filelist)
varnames = fill("",n)
methods = ["Lin24","Log16","R16+Blosc","R14+Blosc","R12+Blosc",
            "R16+LZ4HC","R14+LZ4HC","R12+LZ4HC",
            "ZFP15","ZFP12","ZFP9"]
stats = ["cf","meanerror","absmedian","abs90%","absmax","decmedian","dec90%","decmax"]

E = fill(0f0,n,length(methods),length(stats))

function statscalc(X0::Array{T,N},Xa::Array{T,N}) where {T,N}
    Xamean = mean(Xa)
    X0mean = mean(X0)
    meanerror = abs(X0mean-Xamean)/X0mean       # normalised

    abserr = abs.(X0-Xa)/X0mean                 # normalised
    abserr = sort(vec(abserr))
    abs50 = T(quantile(abserr,.5,sorted=true))
    abs90 = T(quantile(abserr,.9,sorted=true))
    absmax = abserr[end]

    # exception cases for decimal error
    Xa0 = Xa .== zero(T)
    X00 = X0 .== zero(T)

    Xas = signbit.(Xa)
    X0s = signbit.(X0)

    decerr = abs.(log10.(abs.(X0./Xa)))

    decerr[Xas .⊻ X0s] .= Inf        # signs don't match
    decerr[Xa0 .⊻ X00] .= Inf        # one is zero the other one not
    decerr[Xa0 .& X00] .= zero(T)    # both zero, set decimal error to zero

    decerr = sort(vec(decerr))
    dec50 = T(quantile(decerr,.5,sorted=true))
    dec90 = T(quantile(decerr,.9,sorted=true))
    decmax = decerr[end]

    return [meanerror,abs50,abs90,absmax,dec50,dec90,decmax]
end

for (i,file) in enumerate(filelist)
    varname = split(split(file,"_")[end],".")[1]
    varnames[i] = varname
    println("---")
    println("Reading $varname")
    X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    # transpose array
    X = permutedims(X,[3,2,1])  # unravels as longitude, latitude, level

    negatives = any(X .< 0f0)  # check for negative entries

    ## LinQuant24
    E[i,1,1] = 32/24     # compression factor
    Xc = Array{Float32}(LinQuant24Array(X))
    E[i,1,2:end] = statscalc(X,Xc)
    println("LinQuant24 finished.")

<<<<<<< HEAD
    ## LogQuant16
=======
    ## LogQuant16
>>>>>>> 5baed54df0982b9ab5d73d7af6a9f125ceb1e66c
    E[i,2,1] = 32/16
    if negatives       # don't compress
        E[i,2,2:end] .= 0
    else
        Xc = Array{Float32}(LogQuant16Array(X))
        E[i,2,2:end] = statscalc(X,Xc)
    end
    println("LogQuant16 finished.")

<<<<<<< HEAD
    ## Round+lossless
=======
    ## Round+lossless
>>>>>>> 5baed54df0982b9ab5d73d7af6a9f125ceb1e66c
    for (j,r) in enumerate([7,5,3])
        Xc = round(X,r)         # RN16,14,12
        E[i,2+j,2:end] = statscalc(X,Xc)
        E[i,5+j,2:end] = E[i,2+j,2:end]

        Blosc.set_compressor("blosclz")
        E[i,2+j,1] = sizeof(X)/sizeof(compress(Xc,level=5))
        Blosc.set_compressor("lz4hc")
        E[i,5+j,1] = sizeof(X)/sizeof(compress(Xc,level=5))
        println("RN$(r+9)+Blosc/LZ4HC finished.")
    end

<<<<<<< HEAD
    ## Zfp
=======
    ## Zfp
>>>>>>> 5baed54df0982b9ab5d73d7af6a9f125ceb1e66c
    for (j,r) in enumerate([15,12,9])
        Xc = zfp_compress(X,precision=r)
        E[i,8+j,1] = sizeof(X)/sizeof(Xc)
        X2 = zfp_decompress(Xc)
        E[i,8+j,2:end] = statscalc(X,X2)
        println("ZFP$r finished.")
    end

    @save "/Users/milan/cams/error/linlogroundzfp_all.jld" varnames methods stats E
end
