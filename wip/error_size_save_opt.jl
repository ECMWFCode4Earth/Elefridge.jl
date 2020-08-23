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
methods = ["Blosc_L5","Blosc_L9","LZ4HC_L5","LZ4HC_L9","ZFP","ZFP+1","ZFP-1"]
stats = ["cf","meanerror","absmedian","abs90%","absmax","decmedian","dec90%","decmax"]

EE = load("/Users/milan/cams/error/roundzfp_all_opt_gridded.jld")
# E = fill(0f0,n,length(methods),length(stats))
E = EE["E"]

D = load("/Users/milan/cams/entropy/keepbits_98.jld")

## RUN
function rs_zfp(r::Int)
    if r > 6
        rd = 6
    elseif r > 2
        rd = 5
    else
        rd = 4
    end
    return r+rd
end

function statscalc(X0::Array{T,N},Xa::Array{T,N}) where {T,N}
    Xamean = mean(Xa)
    X0mean = mean(X0)
    meanerror = abs(X0mean-Xamean)/X0mean       # normalised

    X0amean = mean(abs.(X0))

    abserr = abs.(X0-Xa)/X0amean                # normalised
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
    varname = split(split(file,"0_")[end],".")[1]
    varnames[i] = varname
    if varname in ["etadot","w","d","vo"]
        println("---")
        println("Reading $varname")
        X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

        # transpose array
        X = permutedims(X,[3,2,1])  # unravels as longitude, latitude, level

        # number of bits of information
        r = D["infbits"][D["varnames"] .== varname][1]-9
        r = max(r,2)    # don't use less than two

        println("$r real information bits.")

        Xr = round(X,r)

        ## Round+BLOSC
        E[i,1,2:end] = statscalc(X,Xr)
        E[i,2,2:end] = E[i,1,2:end]         # copy across for Blosc L9 + LZ4HC L5/9
        E[i,3,2:end] = E[i,1,2:end]
        E[i,4,2:end] = E[i,1,2:end]

        # Blosc.set_compressor("blosclz")
        # E[i,1,1] = sizeof(X)/sizeof(compress(Xr,level=5))
        # E[i,2,1] = sizeof(X)/sizeof(compress(Xr,level=9))
        # Blosc.set_compressor("lz4hc")
        # E[i,3,1] = sizeof(X)/sizeof(compress(Xr,level=5))
        # E[i,4,1] = sizeof(X)/sizeof(compress(Xr,level=9))
        println("Blosc/LZ4HC finished.")

        ## Zfp
        # transpose array for zfp vertical, lat, lon
        X = permutedims(X,[1,2,3])

        r_zfp = rs_zfp(r)
        Xc = zfp_compress(X,precision=r_zfp)
        E[i,5,1] = sizeof(X)/sizeof(Xc)
        X2 = zfp_decompress(Xc)
        E[i,5,2:end] = statscalc(X,X2)

        Xc = zfp_compress(X,precision=r_zfp+1)
        E[i,6,1] = sizeof(X)/sizeof(Xc)
        X2 = zfp_decompress(Xc)
        E[i,6,2:end] = statscalc(X,X2)

        Xc = zfp_compress(X,precision=r_zfp-1)
        E[i,7,1] = sizeof(X)/sizeof(Xc)
        X2 = zfp_decompress(Xc)
        E[i,7,2:end] = statscalc(X,X2)

        println("ZFP$r finished.")
        println(E[i,:,1])
        @save "/Users/milan/cams/error/roundzfp_all_opt_gridded2.jld" varnames methods stats E
    end
end
