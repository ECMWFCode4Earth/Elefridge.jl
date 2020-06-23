using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Blosc
using Elefridge

path = "/Users/milan/cams"
allvars = ["no2","go3","so2","aermr04","aermr05","aermr06","ch4","co"]
vari = allvars[3]

filelist = filter(x->endswith(x,vari*".grib"),readdir(path))
X = xr.open_dataarray(joinpath(path,filelist[1]),engine="cfgrib").data

# volume and error norms (columns):
# vol, mean, median, 90% and max
# total 5

# compression methods (rows):
# LinQuant24, LogQuant16, Bitgroom16/18/20/24 + Blosc, RoundNearest16/18/20/24 + Blosc
# total 10

# storing array
EV = Array{Float32,2}(undef,10,5)

## LinQuant24
EV[1,1] = 32/24     # compression factor
Xc = Array{Float32}(LinQuant24Array(X))
de = abs.(log10.(X ./ Xc))
EV[1,2] = mean(de)
EV[1,3] = median(de)
EV[1,4] = percentile(vec(de),90)
EV[1,5] = maximum(de)

# LoqQuant16
EV[2,1] = 2
Xc = Array{Float32}(LogQuant16Array(X))
de = abs.(log10.(X ./ Xc))
EV[2,2] = mean(de)
EV[2,3] = median(de)
EV[2,4] = percentile(vec(de),90)
EV[2,5] = maximum(de)

## Bitgroom16 + Blosc
Xc = groom(X,7)
de = abs.(log10.(X ./ Xc))
EV[3,2] = mean(de)
EV[3,3] = median(de)
EV[3,4] = percentile(vec(de),90)
EV[3,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[3,1] = sizeof(Xc)/sizeof(Xcc)

# Bitgroom18 + Blosc
Xc = groom(X,9)
de = abs.(log10.(X ./ Xc))
EV[4,2] = mean(de)
EV[4,3] = median(de)
EV[4,4] = percentile(vec(de),90)
EV[4,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[4,1] = sizeof(Xc)/sizeof(Xcc)

# Bitgroom20 + Blosc
Xc = groom(X,11)
de = abs.(log10.(X ./ Xc))
EV[5,2] = mean(de)
EV[5,3] = median(de)
EV[5,4] = percentile(vec(de),90)
EV[5,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[5,1] = sizeof(Xc)/sizeof(Xcc)

# Bitgroom24 + Blosc
Xc = groom(X,15)
de = abs.(log10.(X ./ Xc))
EV[6,2] = mean(de)
EV[6,3] = median(de)
EV[6,4] = percentile(vec(de),90)
EV[6,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[6,1] = sizeof(Xc)/sizeof(Xcc)

# RoundNearest16 + Blosc
Xc = round(X,7)
de = abs.(log10.(X ./ Xc))
EV[7,2] = mean(de)
EV[7,3] = median(de)
EV[7,4] = percentile(vec(de),90)
EV[7,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[7,1] = sizeof(Xc)/sizeof(Xcc)

# RoundNearest18 + Blosc
Xc = round(X,9)
de = abs.(log10.(X ./ Xc))
EV[8,2] = mean(de)
EV[8,3] = median(de)
EV[8,4] = percentile(vec(de),90)
EV[8,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[8,1] = sizeof(Xc)/sizeof(Xcc)

# RoundNearest20 + Blosc
Xc = round(X,11)
de = abs.(log10.(X ./ Xc))
EV[9,2] = mean(de)
EV[9,3] = median(de)
EV[9,4] = percentile(vec(de),90)
EV[9,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[9,1] = sizeof(Xc)/sizeof(Xcc)

# RoundNearest24 + Blosc
Xc = round(X,15)
de = abs.(log10.(X ./ Xc))
EV[10,2] = mean(de)
EV[10,3] = median(de)
EV[10,4] = percentile(vec(de),90)
EV[10,5] = maximum(de)
Xcc = compress(Xc,level=9)
EV[10,1] = sizeof(Xc)/sizeof(Xcc)

# just blosc
Xcc = compress(X,level=9)
bloscc = sizeof(X)/sizeof(Xcc)

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

names = ["LinQuant24","LogQuant16","Bitgroom16+Blosc","Bitgroom18+Blosc",
            "Bitgroom20+Blosc","Bitgroom24+Blosc","RoundNearest16+Blosc",
            "RoundNearest18+Blosc","RoundNearest20+Blosc","RoundNearest24+Blosc"]

markers = ["s","s","d","d","d","d","o","o","o","o"]

fig,ax = subplots(1,1,figsize=(8,6))
ax.set_xscale("log")

ax.plot(theoerror,theosize,"k")
ax.plot(eps(1f0)/2e0,1,"w",markersize=15,marker="*",lw=0,markeredgewidth=1.4,
            markeredgecolor="k",alpha=1,label="Float32",zorder=10)
ax.plot(eps(1f0)/2e0,bloscc,"w",markersize=10,marker="v",lw=0,markeredgewidth=1.4,
            markeredgecolor="k",alpha=0.7,label="Float32+Blosc",zorder=10)

for i in [1,2,7,8,9,10,3,4,5,6]
    # ax.plot(EV[i,[3,5]],[EV[i,1],EV[i,1]],"k",lw=5.5,zorder=-1)
    ax.plot(EV[i,3:4],[EV[i,1],EV[i,1]],colours[i],
            marker=markers[i],markevery=2,
            markeredgecolor="k",markersize=10,
            lw=5,label=names[i])
    ax.plot(EV[i,4:5],[EV[i,1],EV[i,1]],colours[i],lw=2)
end

ax.set_xlim(eps(1f0)/3,1)
ax.set_ylim(0.6,3.5)
ax.set_xlabel("decimal error")
ax.set_ylabel("compression factor")
ax.legend(loc=4,ncol=3,fontsize=9)
ax.set_title(L"NO$_2$ compression",loc="left")

tight_layout()
