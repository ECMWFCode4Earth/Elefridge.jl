using Elefridge
using PyPlot
using PyCall
xr = pyimport("xarray")
using TranscodingStreams, CodecZstd
using JLD

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
methods = ["ZstdL3","ZstdL3+bittranspose",
            "ZstdL10","ZstdL10+bittranspose",
            "ZstdL22","ZstdL22+bittranspose"]
sbits = Array(1:23)
C = fill(0.0,n,length(methods),length(sbits))

for (i,file) in enumerate(filelist)
    varname = split(split(file,"0_")[end],".")[1]
    varnames[i] = varname

    println("---")
    println("Reading $varname")
    X = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    # transpose array
    X = permutedims(X,[3,2,1])  # unravels as longitude, latitude, level
    Xs = X[1:225,1:225,:]       # only 1/8th of the globe

    for (is,s) in enumerate(sbits)
        print(s)
        Xsr = round(Xs,s)
        Xsrt = bittranspose(Xsr)

        Xsr8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xsr)), sizeof(Xsr))
        Xsrt8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xsrt)), sizeof(Xsrt))

        C[i,1,is] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL3,Xsr8))
        C[i,2,is] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL3,Xsrt8))
        println("ZstdL3 done.")
        C[i,3,is] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL10,Xsr8))
        C[i,4,is] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL10,Xsrt8))
        println("ZstdL10 done.")
        C[i,5,is] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL22,Xsr8))
        C[i,6,is] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL22,Xsrt8))
        println("ZstdL22 done.")
    end
    @save "/Users/milan/cams/entropy/compression_transpose_all_zstd_1_8.jld" C varnames methods sbits
end
