using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Blosc
using Elefridge
using TranscodingStreams, CodecZlib, CodecZstd, CodecLz4
using JLD

DeflateCompressorL9 = DeflateCompressor(level=9,windowbits=15)
ZstdCompressorL22 = ZstdCompressor(level=22)

TranscodingStreams.initialize(DeflateCompressorL9)
TranscodingStreams.initialize(ZstdCompressorL22)

path = "/Users/milan/cams"
#allvars = ["no2","go3","so2","aermr04","aermr05","aermr06","ch4","co"]
allvars = ["co2","no","clwc","ciwc","t","w"]

for vari in [allvars[1]]
    println("--------")
    println("Compressing $vari")

    filelist = filter(x->endswith(x,vari*".grib"),readdir(path))
    X = xr.open_dataarray(joinpath(path,filelist[1]),engine="cfgrib").data

    # volume and error norms (columns):
    # vol, mean, median, 90% and max  of decimal error
    # mean, median, 90% and max of abs error
    # total 9

    # compression methods (rows):
    # LinQuant24, LogQuant16, RoundNearest16/24 + Blosc/Deflate/Zstd/LZ4
    # total 10

    # storing array
    EV = Array{Float32,2}(undef,10,9)

    # transpose array
    X = copy(transpose(X))

    ## LinQuant24
    EV[1,1] = 32/24     # compression factor
    Xc = Array{Float32}(LinQuant24Array(X))
    err = abs.(log10.(X ./ Xc))     # decimal error
    EV[1,2] = mean(err)
    EV[1,3] = median(err)
    EV[1,4] = percentile(vec(err),90)
    EV[1,5] = maximum(err)

    err = abs.(X.-Xc)   # absolute error
    EV[1,6] = mean(err)
    EV[1,7] = median(err)
    EV[1,8] = percentile(vec(err),90)
    EV[1,9] = maximum(err)

    # LoqQuant16
    EV[2,1] = 2
    Xc = Array{Float32}(LogQuant16Array(X))
    err = abs.(log10.(X ./ Xc))
    EV[2,2] = mean(err)
    EV[2,3] = median(err)
    EV[2,4] = percentile(vec(err),90)
    EV[2,5] = maximum(err)

    err = abs.(X.-Xc)
    EV[2,6] = mean(err)
    EV[2,7] = median(err)
    EV[2,8] = percentile(vec(err),90)
    EV[2,9] = maximum(err)

    # RoundNearest16 + Blosc
    Xc = round(X,7)
    err = abs.(log10.(X ./ Xc))
    EV[3,2] = mean(err)
    EV[3,3] = median(err)
    EV[3,4] = percentile(vec(err),90)
    EV[3,5] = maximum(err)
    Blosc.set_compressor("blosclz")
    Xcc = compress(Xc,level=9)
    EV[3,1] = sizeof(Xc)/sizeof(Xcc)
    println("RN16+Blosc: CF=$(EV[3,1])")

    err = abs.(X.-Xc)
    EV[3,6] = mean(err)
    EV[3,7] = median(err)
    EV[3,8] = percentile(vec(err),90)
    EV[3,9] = maximum(err)

    # create a UInt8 view
    Xc8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xc)), sizeof(Xc))

    # RoundNearest16 + LZ4HC
    EV[4,2:9] = EV[3,2:9]
    Blosc.set_compressor("lz4hc")
    Xcc = compress(Xc,level=9)
    EV[4,1] = sizeof(Xc)/sizeof(Xcc)
    println("RN16+LZ4HC: CF=$(EV[4,1])")

    # RoundNearest16 + Zstd
    EV[5,2:9] = EV[3,2:9]
    Xcc = transcode(ZstdCompressorL22,Xc8)
    EV[5,1] = sizeof(Xc)/sizeof(Xcc)
    println("RN16+Zstd: CF=$(EV[5,1])")

    # RoundNearest16 + Deflate
    EV[6,2:9] = EV[3,2:9]
    Xcc = transcode(DeflateCompressorL9,Xc8)
    EV[6,1] = sizeof(Xc)/sizeof(Xcc)
    println("RN16+Deflate: CF=$(EV[6,1])")

    # RoundNearest24 + Blosc
    Xc = round(X,15)
    err = abs.(log10.(X ./ Xc))
    EV[7,2] = mean(err)
    EV[7,3] = median(err)
    EV[7,4] = percentile(vec(err),90)
    EV[7,5] = maximum(err)
    Blosc.set_compressor("blosclz")
    Xcc = compress(Xc,level=9)
    EV[7,1] = sizeof(Xc)/sizeof(Xcc)

    err = abs.(X.-Xc)
    EV[7,6] = mean(err)
    EV[7,7] = median(err)
    EV[7,8] = percentile(vec(err),90)
    EV[7,9] = maximum(err)

    # create a UInt8 view
    Xc8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xc)), sizeof(Xc))

    # RoundNearest16 + LZ4
    EV[8,2:9] = EV[7,2:9]
    Blosc.set_compressor("lz4hc")
    Xcc = compress(Xc,level=9)
    EV[8,1] = sizeof(Xc)/sizeof(Xcc)

    # RoundNearest16 + Zstd
    EV[9,2:9] = EV[7,2:9]
    Xcc = transcode(ZstdCompressorL22,Xc8)
    EV[9,1] = sizeof(Xc)/sizeof(Xcc)

    # RoundNearest16 + Deflate
    EV[10,2:9] = EV[7,2:9]
    Xcc = transcode(DeflateCompressorL9,Xc8)
    EV[10,1] = sizeof(Xc)/sizeof(Xcc)

    # just lz4hc
    Blosc.set_compressor("lz4hc")
    Xcc = compress(X,level=9)
    lz4hcc = sizeof(X)/sizeof(Xcc)

    save("cams/error/$vari.jld","EV",EV,"lz4",lz4hcc)
    println("$vari stored.")
end
