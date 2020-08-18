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

ZstdCompressorL5 = ZstdCompressor(level=5)
TranscodingStreams.initialize(ZstdCompressorL5)

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,"_go3.grib"),readdir(path))
grib = xr.open_dataarray(joinpath(path,filelist[end]),engine="cfgrib")
X = grib.data
lat = grib.latitude.data
lon = grib.longitude.data

level = 85
o3 = X[level,:,:]
o3 = copy(o3')      # transpose to have longitude first (better compression)

## compression round + LZ4HC
rbits_ll = [23,7,5,3,1,0]
cfs_ll = fill(0.0,length(rbits_ll))     # compression factors
decerr_ll = fill(0.0,length(rbits_ll))
Blosc.set_compressor("lz4hc")

for (i,r) in enumerate(rbits_ll)
    o3r = bittranspose(round(X,r))
    o3r8 = unsafe_wrap(Array, Ptr{UInt8}(pointer(o3r)), sizeof(o3r))
    o3c = transcode(ZstdCompressorL5,o3r8)
    # o3c = compress(bittranspose(o3r))
    cfs_ll[i] = sizeof(X)/sizeof(o3c)
    # decerr_ll[i] = maximum(vec(abs.(log2.(abs.(X./o3r)))))
end

## compression zfp
rbits_zfp = [23,13,10,8,5,4]               # corresponds to the sigbits from above
cfs_zfp = fill(0.0,length(rbits_zfp))     # compression factors
decerr_zfp = fill(0.0,length(rbits_zfp))   # decimal error

# zfp lossless
cfs_zfp[1] = sizeof(X)/sizeof(zfp_compress(X))

for (i,r) in enumerate(rbits_zfp[2:end])
    local o3c = zfp_compress(X,precision=r)
    cfs_zfp[i+1] = sizeof(X)/sizeof(o3c)
    local o3r = zfp_decompress(o3c)
    decerr_zfp[i+1] = median(vec(abs.(log2.(abs.(X./o3r)))))
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
    o3c = zfp_decompress(zfp_compress(o3,precision=rbits_zfp[i-1]))
    ax1.pcolormesh(lon[lonwhere],lat,o3c[lonwhere,:]',
        vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree(),cmap=cmap)
end

ax0.coastlines(color="k",lw=0.5)
ax1.coastlines(color="k",lw=0.5)

for longi in lond[2:end-1]
    ax0.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
    ax1.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
end

ax0.text(0.02,0.02,"Compression\n               factor",transform=ax0.transAxes,fontsize=8)
ax0.text(0.02,0.92,"         Significant\nbits retained",transform=ax0.transAxes,fontsize=8)

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

tight_layout(rect=[0.0,0.08,1,1])
savefig("/Users/milan/git/Elefridge.jl/maps/zfp_lossless_o3_$level.png",dpi=300)
close(fig)
