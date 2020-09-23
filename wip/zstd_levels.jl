using Elefridge
using ZfpCompression
using PyPlot
using PyCall
using Blosc
using ColorSchemes
using StatsBase, Statistics
using TranscodingStreams, CodecZstd
using BenchmarkTools
xr = pyimport("xarray")
ccrs = pyimport("cartopy.crs")

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,"_t.grib"),readdir(path))
grib = xr.open_dataarray(joinpath(path,filelist[end]),engine="cfgrib")
X = grib.data
lat = grib.latitude.data
lon = grib.longitude.data

X = permutedims(X,[3,2,1])      # for lossless longitude first

## DEFINE COMPRESSORS
ZstdCompressors = Array{ZstdCompressor,1}(undef,22)

for i in 1:22
    ZstdCompressors[i] = ZstdCompressor(level=i)
    TranscodingStreams.initialize(ZstdCompressors[i])
end

## Compress X with round+Zstd
Xr = round(X[1:450,:,:],7)     # temperature contains 7 informationful signficant bits
Xr8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xr)), sizeof(Xr))

cfs = fill(0.0,22)    # compression factors
speed_comp = fill(0.0,22)         # compression speeds and decompression
speed_decomp = fill(0.0,22)

for i in 1:22       # do the first one twice for pre-compile
    println(i)
    Xc = transcode(ZstdCompressors[i],Xr8)
    cfs[i] = sizeof(Xr)/sizeof(Xc)*2/(900*451)*(348528)    # relative to Float64 and unstructured grid
    speed_comp[i] = (@benchmark Xc = transcode(ZstdCompressors[$i],Xr8)).times[1]/1e9
    speed_decomp[i] = (@benchmark Xuc = transcode(ZstdDecompressor,$Xc)).times[1]/1e9
    # @assert Xuc == Xr8      # check whether compression is indeed reversible
end

## Compress with Zfp
Xs = permutedims(X[1:450,:,:],[3,2,1])      # for zfp vert x lat x lon

Xzfp = zfp_compress(Xs,precision=13)
speed_comp_zfp = (@benchmark Xzfp = zfp_compress(Xs,precision=13)).times[1]/1e9
cfs_zfp = sizeof(Xs)/sizeof(Xzfp)*2/(900*451)*(348528)
speed_decomp_zfp = (@benchmark zfp_decompress(Xzfp)).times[1]/1e9


## Compress with LinQuant24
Xlin24 = LinQuant16Array(Xs)
Xlog16 = LogQuant16Array(Xs)

speed_comp_lin24 = (@benchmark Xlin24 = LinQuant16Array(Xs)).times[1]/1e9       # use 16 as UInt24 is slow
speed_comp_log16 = (@benchmark Xlog16 = LogQuant16Array(Xs)).times[1]/1e9

speed_decomp_lin24 = (@benchmark Array{Float32}(Xlin24)).times[1]/1e9
speed_decomp_log16 = (@benchmark Array{Float32}(Xlog16)).times[1]/1e9

## Plotting
ioff()
fig,(ax1,ax2,ax3) = subplots(1,3,sharex=true,figsize=(10,4))

ax1.bar(1:22,cfs,alpha=0.7,color="C1")
ax1.plot([0,23],ones(2)*cfs_zfp,"C2",lw=2)
ax1.plot([0,23],ones(2)*64/24,"grey",lw=2,ls="--")
ax1.plot([0,23],ones(2)*4,"k",lw=2,ls="--")

# speed in MB/s
ax2.plot([0,23],ones(2)*sizeof(Xr)/speed_comp_lin24/1e6,ls="--",lw=2,color="grey",label="LinQuant24")
ax2.plot([0,23],ones(2)*sizeof(Xr)/speed_comp_log16/1e6,ls="--",lw=2,color="k",label="LogQuant16")
ax2.plot([0,23],ones(2)*sizeof(Xr)/speed_comp_zfp/1e6,lw=2,color="C2",label="Zfp")
ax2.bar(1:22,sizeof(Xr)./speed_comp/1e6,alpha=0.7,color="C1",label="Zstd")

ax3.plot([0,23],ones(2)*sizeof(Xr)/speed_decomp_lin24/1e6,ls="--",lw=2,color="grey")
ax3.plot([0,23],ones(2)*sizeof(Xr)/speed_decomp_log16/1e6,ls="--",lw=2,color="k")
ax3.plot([0,23],ones(2)*sizeof(Xr)/speed_decomp_zfp/1e6,lw=2,color="C2")
ax3.bar(1:22,sizeof(Xr)./speed_decomp/1e6,alpha=0.7,color="C1")

ax2.set_yscale("log")
ax3.set_yscale("log")

ax1.set_ylim(1,60)
ax2.set_ylim(1,5000)
ax3.set_ylim(1,5000)
ax1.set_xlim(0,23)

ax1.set_xticks([1,5,10,15,20])
ax1.set_xticks(1:22,minor=true)

ax1.set_yticks(10:10:60)
ax1.set_yticklabels(["10x","20x","30x","40x","50x","60x"])

ax2.set_yticks([1,10,100,1000,5000])
ax2.set_yticklabels([1,10,100,1000,"[MB/s]"])
ax3.set_yticks([1,10,100,1000,5000])
ax3.set_yticklabels([1,10,100,1000,"[MB/s]"])

ax1.set_xlabel("Compression level")
ax2.set_xlabel("Compression level")
ax3.set_xlabel("Compression level")

# ax2.set_ylabel("MB/s")
# ax3.set_ylabel("MB/s")

ax1.set_title("Compression factors",loc="left")
ax2.set_title("Compression speed",loc="left")
ax3.set_title("Decompression speed",loc="left")

ax1.set_title("a",loc="right",fontweight="bold")
ax2.set_title("b",loc="right",fontweight="bold")
ax3.set_title("c",loc="right",fontweight="bold")

ax2.legend(loc=3,framealpha=1)

tight_layout()
savefig("/Users/milan/git/Elefridge.jl/plots/level_speed.png",dpi=100)
close(fig)
