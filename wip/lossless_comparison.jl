using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Blosc
using Elefridge
using TranscodingStreams, CodecZlib, CodecZstd, CodecLz4
using PyPlot

DeflateCompressorL9 = DeflateCompressor(level=9,windowbits=15)
ZstdCompressorL22 = ZstdCompressor(level=22)

TranscodingStreams.initialize(DeflateCompressorL9)
TranscodingStreams.initialize(ZstdCompressorL22)

path = "/Users/milan/cams"
allvars = ["no2","go3","so2","aermr04","aermr05","aermr06","ch4","co"]
vari = allvars[1]

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
Xcc = compress(Xc,level=5)
EV[3,1] = sizeof(Xc)/sizeof(Xcc)

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

# RoundNearest16 + Zstd
EV[5,2:9] = EV[3,2:9]
Xcc = transcode(ZstdCompressorL22,Xc8)
EV[5,1] = sizeof(Xc)/sizeof(Xcc)

# RoundNearest16 + Deflate
EV[6,2:9] = EV[3,2:9]
Xcc = transcode(DeflateCompressorL9,Xc8)
EV[6,1] = sizeof(Xc)/sizeof(Xcc)

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

# just blosc
Blosc.set_compressor("lz4hc")
Xcc = compress(X,level=9)
lz4hcc = sizeof(X)/sizeof(Xcc)

## theoretical
theosize = [32/i for i in 32:-1:5]
x0 = 1f0
theoerror = [log10(nextfloat(x0,2^i)) for i in 0:27]

## plotting
pygui(true)
n = size(EV)[1]

cm1 = ColorMap("viridis")
cm2 = ColorMap("inferno")

function RGB(cv::Tuple)
    r = UInt8(round(cv[1]*255))
    g = UInt8(round(cv[2]*255))
    b = UInt8(round(cv[3]*255))
    return "#"*repr(r)[3:4]*repr(g)[3:4]*repr(b)[3:4]
end

colours = ["0.5","0.3",RGB(cm2(0.9)),RGB(cm2(0.8)),RGB(cm2(0.6)),RGB(cm2(0.4)),
                    RGB(cm1(0.8)),RGB(cm1(0.6)),RGB(cm1(0.4)),RGB(cm1(0.2))]

names = ["LinQuant24","LogQuant16","RoundNearest16+Blosc","RoundNearest16+LZ4HC",
            "RoundNearest16+Zstd","RoundNearest16+Deflate","RoundNearest24+Blosc",
            "RoundNearest24+LZ4HC","RoundNearest24+Zstd","RoundNearest24+Deflate"]

markers = ["s","s","d","d","d","d","o","o","o","o"]

fig,ax = subplots(1,1,figsize=(8,6))
ax.set_xscale("log")

ax.plot(theoerror,theosize,"k")
ax.plot(eps(1f0)/2e0,1,"w",markersize=15,marker="*",lw=0,markeredgewidth=1.4,
            markeredgecolor="k",alpha=1,label="Float32",zorder=10)
ax.plot(eps(1f0)/2e0,lz4hcc,"w",markersize=10,marker="v",lw=0,markeredgewidth=1.4,
            markeredgecolor="k",alpha=0.7,label="Float32+LZ4HC",zorder=10)

for i in [1,2,3,4,5,6,7,8,9,10]
    # ax.plot(EV[i,[3,5]],[EV[i,1],EV[i,1]],"k",lw=5.5,zorder=-1)
    ax.plot(EV[i,3:4],[EV[i,1],EV[i,1]],colours[i],
            marker=markers[i],markevery=2,
            markeredgecolor="k",markersize=10,
            lw=5,label=names[i])
    ax.plot(EV[i,4:5],[EV[i,1],EV[i,1]],colours[i],lw=2)
end

ax.set_xlim(eps(1f0)/3,1)
ax.set_ylim(0.9,4.5)
ax.set_xlabel("decimal error")
ax.set_ylabel("compression factor")
ax.legend(loc=1,ncol=1,fontsize=9)
ax.set_title(L"NO$_2$ compression",loc="left")

ax.text(EV[1,3],EV[1,1],"median   ",rotation=90,ha="center",va="top",fontsize=9)
ax.text(EV[1,4],EV[1,1],"90%   ",rotation=90,ha="center",va="top",fontsize=9)
ax.text(EV[1,5],EV[1,1],"max   ",rotation=90,ha="center",va="top",fontsize=9)

tight_layout()
