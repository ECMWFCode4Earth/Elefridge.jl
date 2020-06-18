using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
using Blosc
using Elefridge

path = "/Users/milan/cams"
allvars = ["no2","go3","so2","aermr04","aermr05","aermr06","ch4","co"]
vari = allvars[1]

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
theosize = [32/i for i in 32:-1:10]
x0 = 1f0
theoerror = [log10(nextfloat(x0,2^i)) for i in 0:22]

## plotting
pygui(true)
n = size(EV)[1]
colours = ["C$i" for i in 0:n-1]
names = ["LinQuant24","LogQuant16","Bitgroom16+Blosc","Bitgroom18+Blosc",
            "Bitgroom20+Blosc","Bitgroom24+Blosc","RoundNearest16+Blosc",
            "RoundNearest18+Blosc","RoundNearest20+Blosc","RoundNearest24+Blosc"]

fig,ax = subplots(1,1,figsize=(8,6))
ax.set_xscale("log")

for i in 1:n
    ax.plot(EV[i,3:4],[EV[i,1],EV[i,1]],colours[i],lw=3,zorder=-1,label=names[i])
    ax.plot(EV[i,4:5],[EV[i,1],EV[i,1]],colours[i],lw=1.5,zorder=-1)
end

ax.scatter(EV[:,2],EV[:,1],50,colours,marker="d",edgecolor="k")
ax.scatter(EV[:,3],EV[:,1],50,colours,marker="o",edgecolor="k")

ax.scatter(eps(1f0)/2e0,1,140,"#00FF00",marker="*",
            edgecolor="k",alpha=0.7,label="Float32",zorder=10)
ax.scatter(eps(1f0)/2e0,bloscc,140,"#00FF00",marker="s",
            edgecolor="k",alpha=0.7,label="Float32+Blosc",zorder=10)
ax.plot(theoerror,theosize,"k")


ax.set_xlim(eps(1f0)/3,1)
ax.set_ylim(0.9,3.5)
ax.set_xlabel("decimal error")
ax.set_ylabel("compression factor")
ax.legend(loc=1)

tight_layout()
