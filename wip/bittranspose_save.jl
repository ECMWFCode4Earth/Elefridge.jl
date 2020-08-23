using Elefridge
using PyPlot
using PyCall
using Blosc
xr = pyimport("xarray")

using TranscodingStreams, CodecZstd

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,"_go3.grib"),readdir(path))
grib = xr.open_dataarray(joinpath(path,filelist[end]),engine="cfgrib")
X = grib.data
lat = grib.latitude.data
lon = grib.longitude.data

X = permutedims(X,[3,2,1])      # for lossless longitude first
Xs = X[1:225,1:225,:]

sbits = Array(1:23)
cfs = fill(0.0,length(sbits))     # compression factors for un-transposed
cft = fill(0.0,length(sbits))     # compression factors for transposed

Blosc.set_compressor("lz4hc")

ZstdCompressorL22 = ZstdCompressor(level=22)
TranscodingStreams.initialize(ZstdCompressorL22)

for (i,s) in enumerate(sbits)
    print(s)
    Xsr = round(Xs,s)
    Xsrt = bittranspose(Xsr)

    # Xsr8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xsr)), sizeof(Xsr))
    # Xsrt8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xsrt)), sizeof(Xsrt))

    # cfs[i] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL22,Xsr8))
    # cft[i] = sizeof(Xs)/sizeof(transcode(ZstdCompressorL22,Xsrt8))

    cfs[i] = sizeof(Xs)/sizeof(compress(Xsr,level=9))
    cft[i] = sizeof(Xs)/sizeof(compress(Xsrt,level=9))

end

## PLOT

plot(cfs./cft)
xlabel("significant bits")
ylabel("compression factor: untransposed/transposed")
