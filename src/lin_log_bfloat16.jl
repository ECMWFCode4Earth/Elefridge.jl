using NetCDF, LinearAlgebra
using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
#pygui(:qt)
using PyPlot
using MultivariateStats
# eofs = pyimport("eofs")
using BFloat16s
using Sonums

path = "/Users/milan/Downloads/cams/ieee"
filelist = filter(x->endswith(x,"no2.grib"),readdir(path))
gribfile = xr.open_dataset(joinpath(path,filelist[1]),engine="cfgrib")
no2 = gribfile.no2.values
level = 1:size(no2)[1]

# train sonums
trainSonum16(no2)

no2lin24 = Array{Float32}(LinQuant24Array(no2))
no2log16 = Array{Float32}(LogQuant16Array(no2))
no2bf16 = Float32.(BFloat16.(no2))
no2s16 = Float32.(Sonum16.(no2))

## error 24bit
de = abs.(log10.(no2 ./ no2lin24))

de24m = mean(de,dims=2)[:,1]
de24mi = median(de,dims=2)[:,1]
de24min = minimum(de,dims=2)[:,1]
de24max = maximum(de,dims=2)[:,1]
de24p10 = [percentile(vec(de[i,:]),10) for i in level]
de24p90 = [percentile(vec(de[i,:]),90) for i in level]

# and with 16
de = abs.(log10.(no2 ./ no2log16))

de16m = mean(de,dims=2)[:,1]
de16mi = median(de,dims=2)[:,1]
de16min = minimum(de,dims=2)[:,1]
de16max = maximum(de,dims=2)[:,1]
de16p10 = [percentile(vec(de[i,:]),10) for i in level]
de16p90 = [percentile(vec(de[i,:]),90) for i in level]

# and for bfoat16 = bitgrooming
de = abs.(log10.(no2 ./ no2bf16))

debf16m = mean(de,dims=2)[:,1]
debf16mi = median(de,dims=2)[:,1]
debf16min = minimum(de,dims=2)[:,1]
debf16max = maximum(de,dims=2)[:,1]
debf16p10 = [percentile(vec(de[i,:]),10) for i in level]
debf16p90 = [percentile(vec(de[i,:]),90) for i in level]

# and for bfoat16 = bitgrooming
de = abs.(log10.(no2 ./ no2s16))

des16m = mean(de,dims=2)[:,1]
des16mi = median(de,dims=2)[:,1]
des16min = minimum(de,dims=2)[:,1]
des16max = maximum(de,dims=2)[:,1]
des16p10 = [percentile(vec(de[i,:]),10) for i in level]
des16p90 = [percentile(vec(de[i,:]),90) for i in level]

## plot
fig,ax = subplots(1,1,figsize=(6,10))
ax.invert_yaxis()

ax.semilogx(de24m,level,"C0",lw=3,label="24-bit lin")
ax.semilogx(de16m,level,"C1",lw=3,label="16-bit log")
ax.semilogx(debf16m,level,"C2",lw=3,label="bfloat16")
ax.semilogx(des16m,level,"C3",lw=3,label="15-bit max-entropy")

ax.semilogx(de24p90,level,"C0",lw=2)
ax.semilogx(de16p90,level,"C1",lw=2)
ax.semilogx(debf16p90,level,"C2",lw=2)
ax.semilogx(des16p90,level,"C3",lw=2)

ax.semilogx(de24max,level,"C0",lw=1)
ax.semilogx(de16max,level,"C1",lw=1)
ax.semilogx(debf16max,level,"C2",lw=1)
ax.semilogx(des16max,level,"C3",lw=1)

ax.plot(0,0,"grey",lw=3,label="mean")
ax.plot(0,0,"grey",lw=2,label="90%")
ax.plot(0,0,"grey",lw=1,label="max")

ax.set_title(L"NO$_2$ compression error")
ax.set_ylabel("model level")
ax.set_xlabel("decimal error")
ax.set_ylim(137,1)
ax.legend(loc=3,ncol=2)

# a‚x.set_xlim(1e-7,1)
tight_layout()
