using Elefridge
using ZfpCompression
using PyPlot
using PyCall
using ColorSchemes
using StatsBase, Statistics
xr = pyimport("xarray")
ccrs = pyimport("cartopy.crs")

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,"_go3.grib"),readdir(path))
grib = xr.open_dataarray(joinpath(path,filelist[end]),engine="cfgrib")
X = grib.data
lat = grib.latitude.data
lon = grib.longitude.data

level = 85
X = permutedims(X,[3,2,1])      # for lossless longitude first
o3 = X[:,:,level]

## compression round
o3r = round(o3,5)

## compression zfp
X = permutedims(X,[3,2,1])      # for zfp vert x lat x lon
o3z = zfp_decompress(zfp_compress(X,precision=10))[level,:,:]'

## PLOT
lat_div(n::Integer) = Array(-90:180/(n-1):90)
lon_div(n::Integer) = Array(0:360/n:360)
lond = lon_div(6)

# pygui(true)
# ion()
ioff()
cmap = ColorMap(ColorSchemes.turku.colors)
vmin,vmax = extrema(o3)
fig,(ax0,ax1,ax2) = plt.subplots(3,1,figsize=(7,10),sharex=true,sharey=true,
                    subplot_kw=Dict("projection"=>ccrs.Robinson(central_longitude=180)))

ax2.pcolormesh(lon,lat,o3',transform=ccrs.PlateCarree();vmin,vmax,cmap)
ax1.pcolormesh(lon,lat,o3',transform=ccrs.PlateCarree();vmin,vmax,cmap)
ax0.pcolormesh(lon,lat,o3',transform=ccrs.PlateCarree();vmin,vmax,cmap)

# ax0.coastlines(color="k",lw=0.5)
# ax1.coastlines(color="k",lw=0.5)
#
# for longi in lond[2:end-1]
#     ax0.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
#     ax1.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
# end

ax0.text(0.9,0.95,"a",fontweight="bold",transform=ax0.transAxes)
ax1.text(0.9,0.95,"b",fontweight="bold",transform=ax1.transAxes)
ax2.text(0.9,0.95,"c",fontweight="bold",transform=ax2.transAxes)

tight_layout(rect=[0,0,1,1])
savefig("/Users/milan/git/Elefridge.jl/maps/turing2_o3_$level.png",dpi=200)
close(fig)
