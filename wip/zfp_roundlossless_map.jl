using Elefridge
using ZfpCompression
using PyPlot
using PyCall
using Blosc
using ColorSchemes
using StatsBase, Statistics
using TranscodingStreams, CodecZstd
xr = pyimport("xarray")
ccrs = pyimport("cartopy.crs")

ZstdCompressorL3 = ZstdCompressor(level=3)
TranscodingStreams.initialize(ZstdCompressorL3)

ZstdCompressorL10 = ZstdCompressor(level=10)
TranscodingStreams.initialize(ZstdCompressorL10)

ZstdCompressorL22 = ZstdCompressor(level=22)
TranscodingStreams.initialize(ZstdCompressorL22)

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,"_t.grib"),readdir(path))
grib = xr.open_dataarray(joinpath(path,filelist[end]),engine="cfgrib")
X = grib.data
lat = grib.latitude.data
lon = grib.longitude.data

level = 85
X = permutedims(X,[3,2,1])      # for lossless longitude first
o3 = X[:,:,level]

## compression round + LZ4HC/Zstd
rbits_ll = [23,9,8,7,6,5]
cfs_ll = fill(0.0,length(rbits_ll))     # compression factors
decerr_ll = fill(0.0,length(rbits_ll))
# Blosc.set_compressor("lz4hc")

for (i,r) in enumerate(rbits_ll)
    Xr = round(X[1:450,:,:],r)      # use only half the longitudes for speed-up

    Xr8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(Xr)), sizeof(Xr))
    Xc = transcode(ZstdCompressorL3,Xr8)

    # Xc = compress(Xr,level=5)
    cfs_ll[i] = sizeof(Xr)/sizeof(Xc)*2/(900*451)*(348528)    # relative to Float64 and unstructured grid
    # decerr_ll[i] = maximum(vec(abs.(log2.(abs.(X./o3r)))))
end

## compression zfp
function rs_zfp(r::Int)
    if r > 6
        rd = 6
    elseif r > 2
        rd = 5
    else
        rd = 4
    end
    return r+rd
end

X = grib.data      # for zfp vert x lat x lon
rbits_zfp = rs_zfp.(rbits_ll)               # corresponds to the sigbits from above
cfs_zfp = fill(0.0,length(rbits_zfp))       # compression factors
# decerr_zfp = fill(0.0,length(rbits_zfp))  # decimal error

O3plot = fill(0f0,length(rbits_zfp)-1,size(o3)...)

# zfp lossless, relative to Float64 and unstructured grid
cfs_zfp[1] = sizeof(X)/sizeof(zfp_compress(X))*2/(900*451)*(348528)

for (i,r) in enumerate(rbits_zfp[2:end])
    Xc = zfp_compress(X,precision=r)
    cfs_zfp[i+1] = sizeof(X)/sizeof(Xc)*2/(900*451)*(348528)
    O3plot[i,:,:] = zfp_decompress(Xc)[level,:,:]'
    # local o3r = zfp_decompress(o3c)
    # decerr_zfp[i+1] = median(vec(abs.(log2.(abs.(X./o3r)))))
end

lat_div(n::Integer) = Array(-90:180/(n-1):90)
lon_div(n::Integer) = Array(0:360/n:360)

## PLOT
ioff()
cmap = ColorMap(ColorSchemes.turku.colors)
vmin,vmax = extrema(o3)
fig,(ax0,ax1) = plt.subplots(2,1,figsize=(7,8),sharex=true,sharey=true,
                    subplot_kw=Dict("projection"=>ccrs.Robinson(central_longitude=180)))

pos = ax1.get_position()
cax = fig.add_axes([pos.x0,0.07,pos.x1-pos.x0,0.02])

lond = lon_div(length(rbits_ll))

# Float32 / lossless
scale = 1e7
pcm = ax0.pcolormesh(lon[lon .<= lond[2]],lat,scale*o3[lon .<= lond[2],:]',
    vmin=scale*vmin,vmax=scale*vmax,transform=ccrs.PlateCarree(),cmap=cmap)
ax1.pcolormesh(lon[lon .<= lond[2]],lat,o3[lon .<= lond[2],:]',
    vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree(),cmap=cmap)

cbar = colorbar(pcm,cax=cax,orientation="horizontal")
cbar.set_label(L"mixing ratio [10$^{-7}$ kg/kg]")
cbar.set_ticks([1,2,3])

# round + lossless
for i in 3:length(lond)
    lonwhere = (lon .<= lond[i]) .& (lon .> lond[i-1]-0.5)
    ax0.pcolormesh(lon[lonwhere],lat,round(o3,rbits_ll[i-1])[lonwhere,:]',
        vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree(),cmap=cmap)
end

# zfp
for i in 3:length(lond)
    lonwhere = (lon .<= lond[i]) .& (lon .> lond[i-1]-0.5)
    ax1.pcolormesh(lon[lonwhere],lat,O3plot[i-2,lonwhere,:]',
        vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree(),cmap=cmap)
end

ax0.coastlines(color="k",lw=0.5)
ax1.coastlines(color="k",lw=0.5)

for longi in lond[2:end-1]
    ax0.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
    ax1.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
end

ax0.text(0.02,0.02,"Compression\n               factor",transform=ax0.transAxes,fontsize=8)
ax0.text(0.02,0.93,"         Significant\nbits retained",transform=ax0.transAxes,fontsize=8)
ax0.text(0.88,0.87,"Real information\n       preserved",transform=ax0.transAxes,fontsize=8)

info_preserved = ["100%","99.9%","99%","95%","82%","71%"]

for (i,x) in enumerate([0.16,0.28,0.41,0.54,0.65,0.77])
    ax0.text(x,0.89,info_preserved[i],transform=ax0.transAxes,fontweight="bold",
                fontsize=10,color="white")
    ax1.text(x,0.89,info_preserved[i],transform=ax1.transAxes,fontweight="bold",
                fontsize=10,color="white")
end

for (i,x) in enumerate([0.23,0.33,0.43,0.53,0.63,0.74])
    c = Int(round(cfs_ll[i]))
    r = rbits_ll[i]
    ax0.text(x,0.03,"$(c)x",transform=ax0.transAxes,fontweight="bold",fontsize=12,
                            color="white")
    ax0.text(x,0.98,"$(r)",transform=ax0.transAxes,fontweight="bold",fontsize=10,
                            color="white",va="top")
    ax1.text(x,0.98,"$(r)",transform=ax1.transAxes,fontweight="bold",fontsize=10,
                            color="white",va="top")
    c = Int(round(cfs_zfp[i]))
    ax1.text(x,0.03,"$(c)x",transform=ax1.transAxes,fontweight="bold",fontsize=12,
                            color="white")
end

ax0.set_title(L"O$_3$ compression, round+lossless")
ax0.set_title("a",fontweight="bold",loc="right")

ax1.set_title(L"O$_3$ zfp compression")
ax1.set_title("b",fontweight="bold",loc="right")

tight_layout(rect=[0.02,0.08,1,1])
savefig("/Users/milan/git/Elefridge.jl/maps/zfp_t_$level.png",dpi=200)
close(fig)
