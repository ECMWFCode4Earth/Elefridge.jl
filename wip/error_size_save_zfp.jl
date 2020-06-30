using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Elefridge
using ZfpCompression
using JLD

path = "/Users/milan/cams"
allvars = ["go3","so2","aermr04","aermr05","aermr06","ch4","co","co2","no"]
# allvars = ["no2"]

for vari in allvars
    println("--------")
    println("Compressing $vari")

    filelist = filter(x->endswith(x,vari*".grib"),readdir(path))
    X = xr.open_dataarray(joinpath(path,filelist[1]),engine="cfgrib").data

    # volume and error norms (columns):
    # vol, mean, median, 90% and max  of decimal error
    # mean, median, 90% and max of abs error
    # total 9

    # compression methods (rows):
    # ZFP lossless, ZFP tol = 1e-7,1e-5,1e-11, precision=14,17,20
    # ZFP lossless with RoundNearest16,RoundNearest24
    # total 9

    # storing array
    EV = Array{Float32,2}(undef,9,9)

    ## LOSSLESS
    Xc = zfp_compress(X)
    EV[1,1] = sizeof(X)/sizeof(Xc)     # compression factor
    EV[1,2:end] .= 0

    println("Lossless, $(EV[1,1])")

    # ZFP tol
    for (i,tol) in enumerate([1e-7,1e-9,1e-11])
        Xc = zfp_compress(X,tol=tol)
        EV[1+i,1] = sizeof(X)/sizeof(Xc)     # compression factor

        X2 = similar(X)
        zfp_decompress!(X2,Xc,tol=tol)

        if any(X2 .< 0)
            err = fill(-1,size(X))
        else
            err = abs.(log10.(X ./ X2))         # decimal error
        end
        EV[1+i,2] = mean(err)
        EV[1+i,3] = median(err)
        EV[1+i,4] = percentile(vec(err),90)
        EV[1+i,5] = maximum(err)

        err = abs.(X.-X2)   # absolute error
        EV[1+i,6] = mean(err)
        EV[1+i,7] = median(err)
        EV[1+i,8] = percentile(vec(err),90)
        EV[1+i,9] = maximum(err)
        println("tol=$tol, cf=$(EV[1+i,1]), max dec=$(EV[1+i,5]), max abs=$(EV[1+i,9])")
    end

    # ZFP precision
    for (i,prec) in enumerate([14,17,20])
        Xc = zfp_compress(X,precision=prec)
        EV[4+i,1] = sizeof(X)/sizeof(Xc)     # compression factor

        X2 = similar(X)
        zfp_decompress!(X2,Xc,precision=prec)

        if any(X2 .< 0)
            err = fill(-1,size(X))
        else
            err = abs.(log10.(X ./ X2))         # decimal error
        end
        EV[4+i,2] = mean(err)
        EV[4+i,3] = median(err)
        EV[4+i,4] = percentile(vec(err),90)
        EV[4+i,5] = maximum(err)

        err = abs.(X.-X2)   # absolute error
        EV[4+i,6] = mean(err)
        EV[4+i,7] = median(err)
        EV[4+i,8] = percentile(vec(err),90)
        EV[4+i,9] = maximum(err)
        println("precision=$prec, cf=$(EV[4+i,1]), max dec=$(EV[4+i,5]), max abs=$(EV[4+i,9])")
    end

    # RoundNearest16/24 + zfp lossless
    for (i,bits) in enumerate([7,15])
        X2 = round(X,bits)

        Xc = zfp_compress(X2)
        EV[7+i,1] = sizeof(X)/sizeof(Xc)     # compression factor

        err = abs.(log10.(X ./ X2))         # decimal error
        EV[7+i,2] = mean(err)
        EV[7+i,3] = median(err)
        EV[7+i,4] = percentile(vec(err),90)
        EV[7+i,5] = maximum(err)

        err = abs.(X.-X2)   # absolute error
        EV[7+i,6] = mean(err)
        EV[7+i,7] = median(err)
        EV[7+i,8] = percentile(vec(err),90)
        EV[7+i,9] = maximum(err)
        println("RN=$bits, cf=$(EV[7+i,1]), max dec=$(EV[7+i,5]), max abs=$(EV[7+i,9])")
    end

    save(joinpath(path,"error/$(vari)_zfp.jld"),"EV",EV)
    println("$vari stored.")
end
