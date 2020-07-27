using Elefridge
using ZfpCompression
using PyPlot
using PyCall
xr = pyimport("xarray")
ccrs = pyimport("cartopy.crs")

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,"_go3.grib"),readdir(path))
grib = xr.open_dataarray(joinpath(path,filelist[end]),engine="cfgrib")
X = grib.data
lat = grib.latitude.data
lon = grib.longitude.data

## loop

for level in 5:5:137
    o3 = X[level,:,:]
    o3 = copy(o3')

    ## compression
    rbits = [14,11,9,8,5]
    cfs = fill(0.0,length(rbits))     # compression factors

    for (i,r) in enumerate(rbits)
        o3c = zfp_compress(X,precision=r)
        cfs[i] = sizeof(X)/sizeof(o3c)
    end

    lat_div(n::Integer) = Array(-90:180/(n-1):90)
    lon_div(n::Integer) = Array(0:360/n:360)

    ## PLOT
    ioff()
    vmin,vmax = extrema(o3)
    fig,ax1 = plt.subplots(1,1,figsize=(10,6),sharex=true,sharey=true,
                        subplot_kw=Dict("projection"=>ccrs.Robinson(central_longitude=180)))

    lond = lon_div(length(rbits)+1)

    # Float32
    ax1.pcolormesh(lon[lon .<= lond[2]],lat,o3[lon .<= lond[2],:]',vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree())

    for i in 3:length(lond)
        lonwhere = (lon .<= lond[i]) .& (lon .> lond[i-1]-0.5)
        o3c = zfp_decompress(zfp_compress(o3,precision=rbits[i-2]))
        ax1.pcolormesh(lon[lonwhere],lat,o3c[lonwhere,:]',vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree())
    end

    ax1.coastlines(color="k",lw=0.5)


    for longi in lond[2:end-1]
        ax1.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
        # ax1.plot(repeat([longi,],6),lat_div(31)[end-5:end],"w",transform=ccrs.PlateCarree(),lw=1.5)
    end

    ax1.text(0.127,0.105,"Compression\nfactor",transform=ax1.transAxes,fontweight="bold",color="white")

    for (i,x) in enumerate([0.23,0.33,0.43,0.54,0.64,0.76])
        c = i == 1 ? 1 : Int(round(cfs[i-1]))
        ax1.text(x,0.03,"$(c)x",transform=ax1.transAxes,fontweight="bold",fontsize=14,
                                color="white")

    end

    ax1.set_title(L"O$_3$ zfp compression",fontweight="bold")

    tight_layout()
    savefig("/Users/milan/git/Elefridge.jl/maps/o3/zfp_precision3d_o3_$level.png",dpi=200)
    close(fig)
end
