using Statistics, StatsBase, NetCDF
using JLD2, FileIO
using ZfpCompression, LinLogQuantization, BitInformation
using TranscodingStreams, CodecZstd

ZstdCompressorL3 = ZstdCompressor(level=3)
TranscodingStreams.initialize(ZstdCompressorL3)

ZstdCompressorL10 = ZstdCompressor(level=10)
TranscodingStreams.initialize(ZstdCompressorL10)

ZstdCompressorL22 = ZstdCompressor(level=22)
TranscodingStreams.initialize(ZstdCompressorL22)

path = "/network/aopp/chaos/pred/kloewer/esowc/cams/"
filelist = filter(x->endswith(x,"_v3.nc"),readdir(path))
nvars = length(filelist)

function statscalc( X0::Array{T,N},
                    Xa::Array{T,N},
                    qs=[.5,.67,.8,.9,.95,1]) where {T,N}
    
    # normalised absolute error
    X0mean = mean(X0[X0 .> 0])
    abserr = abs.(X0-Xa)/X0mean
    abserr = vec(abserr)
    sort!(abserr)
    abserrs = [quantile(abserr,q,sorted=true) for q in qs]

    # exception cases for decimal error
    Xa0 = Xa .== zero(T)
    X00 = X0 .== zero(T)

    Xas = signbit.(Xa)
    X0s = signbit.(X0)

    decerr = abs.(log10.(abs.(X0./Xa)))

    decerr[Xas .⊻ X0s] .= Inf        # signs don't match
    decerr[Xa0 .⊻ X00] .= Inf        # one is zero the other one not
    decerr[Xa0 .& X00] .= zero(T)    # both zero, set decimal error to zero

    decerr = vec(decerr)
    sort!(decerr)
    decerrs = [quantile(decerr,q,sorted=true) for q in qs]

    return abserrs,decerrs
end

@load joinpath(path,"../analysis/bits_to_retain.jld2") infbits infbits100

methods = ["Lin16","Log16","RL100","RL99","ZFP100","ZFP99"]
nmethods = length(methods)
nstat = 2   # abs error & dec error
nqs = 6     # number of quantiles

@load joinpath(path,"../analysis","compress_errors.jld2") E C varnames
# varnames = fill("",nvars)
# E = fill(0.0,nvars,nmethods,nstat,nqs)    # error array
# C = fill(0.0,nvars,nmethods)              # compression factor array

for (i,file) in enumerate(filelist)
    varname = split(split(file,"cams_")[end],"_2019")[1]
    varnames[i] = varname
    
    if i > 25
    
        print("$varname, $(infbits[i]-9)bits ")
        ncfile = NetCDF.open(joinpath(path,file))

        # find the variable name by size
        var = [var for var in ncfile.vars if prod(size(var[2])) == 900*451*137][1][1]
        X = ncfile.vars[var][:,:,:]  
        orisize = 2*sizeof(X)     # relative to Float64
        print("L")

        # set negatives to zero (only in aerlg present)
        X[X .< 0] .= 0

        # LinQuant24
        Lin24 = Array{Float32}(LinQuant24Array(X,3))
        abserror,decerror = statscalc(X,Lin24)
        E[i,1,1,:] = abserror
        E[i,1,2,:] = decerror
        C[i,1] = 64/24
        print("a")

        # LogQuant16
        Log16 = Array{Float32}(LogQuant16Array(X,3))
        abserror,decerror = statscalc(X,Log16)
        E[i,2,1,:] = abserror
        E[i,2,2,:] = decerror
        C[i,2] = 64/16
        print("b")

        # Round+lossless 100% preserved information
        Xr = round(X,infbits100[i]-9)   
        abserror,decerror = statscalc(X,Xr)
        E[i,3,1,:] = abserror
        E[i,3,2,:] = decerror
        Xr8 = copy(reinterpret(UInt8,vec(Xr)))
        C[i,3] = orisize/sizeof(transcode(ZstdCompressorL10,Xr8))
        print("c")

        # Round+lossless 99% preserved information
        Xr = round(X,infbits[i]-9)   
        abserror,decerror = statscalc(X,Xr)
        E[i,4,1,:] = abserror
        E[i,4,2,:] = decerror
        Xr8 = copy(reinterpret(UInt8,vec(Xr)))
        C[i,4] = orisize/sizeof(transcode(ZstdCompressorL10,Xr8))
        print("d")

        # ZFP 99% preserved information
        median_abserror = E[i,4,1,1]              # determine precision by matching median abserror
        median_abserror_zfp = median_abserror     # start larger for while loop
        prec = 5                                  # start low and work up
        first_while = true
        print("e")
        while first_while || median_abserror_zfp > median_abserror
            if all(X .> 0)
                logX = log.(X)                        # apply log-preprocessing
                Xc = zfp_compress(logX,precision=prec,nthreads=8)
                C[i,6] = orisize/sizeof(Xc)
                Xz = exp.(zfp_decompress(Xc))
            else
                Xc = zfp_compress(X,precision=prec,nthreads=8)
                C[i,6] = orisize/sizeof(Xc)
                Xz = zfp_decompress(Xc)
            end
            abserror,decerror = statscalc(X,Xz)
            E[i,6,1,:] = abserror
            E[i,6,2,:] = decerror
            median_abserror_zfp = E[i,6,1,1]
            prec += 1
            first_while = false
            print(".")
        end

        # ZFP 100% preserved information
        median_abserror = E[i,3,1,1]              # determine precision by matching median abserror
        first_while = true
        print("f")
        while first_while || median_abserror_zfp > median_abserror
            if all(X .> 0)
                logX = log.(X)                        # apply log-preprocessing
                Xc = zfp_compress(logX,precision=prec,nthreads=8)
                C[i,5] = orisize/sizeof(Xc)
                Xz = exp.(zfp_decompress(Xc))
            else
                Xc = zfp_compress(X,precision=prec,nthreads=8)
                C[i,5] = orisize/sizeof(Xc)
                Xz = zfp_decompress(Xc)
            end
            abserror,decerror = statscalc(X,Xz)
            E[i,5,1,:] = abserror
            E[i,5,2,:] = decerror
            median_abserror_zfp = E[i,5,1,1]
            prec += 1
            first_while = false
            print(".")
        end
        println("; ")
    
        @save joinpath(path,"../analysis","compress_errors.jld2") varnames E C methods
   end
end
                