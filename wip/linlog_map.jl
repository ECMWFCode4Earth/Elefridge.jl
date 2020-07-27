using Elefridge
using ZfpCompression
using PyPlot
using PyCall
using Blosc
xr = pyimport("xarray")
ccrs = pyimport("cartopy.crs")

path = "/Users/milan/cams/gridded/"
filelist = filter(x->endswith(x,"_go3.grib"),readdir(path))
grib = xr.open_dataarray(joinpath(path,filelist[end]),engine="cfgrib")
X = grib.data
lat = grib.latitude.data
lon = grib.longitude.data

o3clin24 = Array(LinQuant24Array(X))
o3clog16 = Array(LogQuant16Array(X))
o3clin16 = Array(LinQuant16Array(X))
o3clog8 = Array(LogQuant8Array(X))
o3clin8 = Array(LinQuant8Array(X))

for level in 10:10:137

    ## compression
    rbits = [10,7,5,3,1]
    cfs = [32/24,32/16,32/16,32/8,32/8]

    lat_div(n::Integer) = Array(-90:180/(n-1):90)
    lon_div(n::Integer) = Array(0:360/n:360)

    ## PLOT
    ioff()
    o3 = X[level,:,:]'
    vmin,vmax = extrema(o3)
    fig,ax1 = plt.subplots(1,1,figsize=(10,6),sharex=true,sharey=true,
                        subplot_kw=Dict("projection"=>ccrs.Robinson(central_longitude=180)))

    lond = lon_div(length(rbits)+1)

    # Float32
    ax1.pcolormesh(lon[lon .<= lond[2]],lat,o3[lon .<= lond[2],:]',vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree())

    for (i,o3c) in zip(3:length(lond),[o3clin24,o3clog16,o3clin16,o3clog8,o3clin8])
        lonwhere = (lon .<= lond[i]) .& (lon .> lond[i-1]-0.5)

        o3cl = o3c[level,:,:]'

        ax1.pcolormesh(lon[lonwhere],lat,o3cl[lonwhere,:]',vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree())
    end

    ax1.coastlines(color="k",lw=0.5)


    for longi in lond[2:end-1]
        ax1.plot(repeat([longi,],31),lat_div(31),"w",transform=ccrs.PlateCarree(),lw=0.25)
        # ax1.plot(repeat([longi,],6),lat_div(31)[end-5:end],"w",transform=ccrs.PlateCarree(),lw=1.5)
    end

    ax1.text(0.127,0.105,"Compression\nfactor",transform=ax1.transAxes,fontweight="bold",color="white")

    for (i,x) in enumerate([0.23,0.33,0.43,0.54,0.64,0.76])
        c = i == 1 ? 1 : Int(round(cfs[i-1]))
        c = i == 2 ? "1.3" : c
        ax1.text(x,0.03,"$(c)x",transform=ax1.transAxes,fontweight="bold",fontsize=14,
                                color="white")
    end

    for (x,txt) in zip([0.24,0.34,0.43,0.53,0.62,0.72],
            ["32-bit\nFloat","24-bit\nLin","16-bit\nLog","16-bit\nLin",
            "8-bit\nLog","8-bit\nLin"])
        ax1.text(x,-0.07,txt,transform=ax1.transAxes)
    end

    ax1.set_title(L"O$_3$ compression: Lin/LogQuantization",fontweight="bold")

    tight_layout()
    savefig("/Users/milan/git/Elefridge.jl/maps/o3/linlog_o3_$level.png",dpi=200)
    close(fig)
end
