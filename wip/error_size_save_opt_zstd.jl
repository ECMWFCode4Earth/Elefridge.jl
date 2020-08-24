using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Blosc
using Elefridge
using JLD
using ZfpCompression
using TranscodingStreams,CodecZstd

ZstdCompressorL3 = ZstdCompressor(level=3)
TranscodingStreams.initialize(ZstdCompressorL3)

ZstdCompressorL10 = ZstdCompressor(level=10)
TranscodingStreams.initialize(ZstdCompressorL10)

ZstdCompressorL22 = ZstdCompressor(level=22)
TranscodingStreams.initialize(ZstdCompressorL22)

path = "/Users/milan/cams/gridded"
filelist = filter(x->endswith(x,".grib"),readdir(path))

n = length(filelist)
varnames = fill("",n)
methods = ["ZstdL3","ZstdL10","ZstdL22"]
stats = ["cf","meanerror","absmedian","abs90%","absmax","decmedian","dec90%","decmax"]

E = load("/Users/milan/cams/error/round_zstd_all_gridded.jld")
E = E["E"]

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

    println("---")
    println("Reading $varname")
    X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    # transpose array
    # use only half the globe
    X = permutedims(X[:,:,1:450],[3,2,1])  # unravels as longitude, latitude, level

    # number of bits of information
    r = D["infbits"][D["varnames"] .== varname][1]-9
    r = max(r,2)    # don't use less than two

    println("$r real information bits.")

    Xr = round(X,r)

    ## Round+Zstd
    E[i,1,2:end] = statscalc(X,Xr)
    E[i,2,2:end] = E[i,1,2:end]         # copy across for Zstd
    E[i,3,2:end] = E[i,1,2:end]

    Xr8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xr)), sizeof(Xr))

    E[i,1,1] = sizeof(X)/sizeof(transcode(ZstdCompressorL3,Xr8))
    println("ZstdL3 finished.")
    E[i,2,1] = sizeof(X)/sizeof(transcode(ZstdCompressorL10,Xr8))
    println("ZstdL10 finished.")
    E[i,3,1] = sizeof(X)/sizeof(transcode(ZstdCompressorL22,Xr8))
    println("ZstdL22 finished.")

    println(E[i,:,1])
    @save "/Users/milan/cams/error/round_zstd_all_gridded.jld" varnames methods stats E
end
