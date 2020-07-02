using NetCDF
using Elefridge
using ZfpCompression
using PyPlot

path = "/Users/milan/cams/nc/"
ncfile = NetCDF.open(joinpath(path,"z_cams_c_ecmf_20200702000000_prod_fc_ml_000_go3.nc"))

no2 = ncfile.vars["go3"][:,:,:,1]
lat = ncfile.vars["latitude"][:]
lon = ncfile.vars["longitude"][:]

# size of array in byte on native octahedral grid for Float32
s = 1*348528*4
level = 100
no2c = zfp_compress(no2[:,:,level],precision=8)
no2d = zfp_decompress(no2c)

maximum(abs.(log10.(no2d./no2[:,:,level])))
maximum(abs.(no2d - no2[:,:,level]))

cf = s/sizeof(no2c)

##
ioff()
fig,(ax1,ax2) = subplots(2,1,sharex=true,sharey=true,figsize=(7,7))

vmin = minimum(no2d)
vmax = maximum(no2d)

ax1.pcolormesh(lon,lat,no2[:,:,level]',vmin=vmin,vmax=vmax)
ax2.pcolormesh(lon,lat,no2d',vmin=vmin,vmax=vmax)

ax1.set_title(L"O$_3$, level "*"$level",loc="left")
ax1.set_title("a",loc="right",fontweight="bold")
ax2.set_title(L"O$_3$, level "*"$level, zfp compression,"*" cf=$(round(cf,digits=1))",loc="left")
ax2.set_title("b",loc="right",fontweight="bold")

ax1.set_xticks([0,90,180,270,360])
ax1.set_yticks([-60,-30,0,30,60])
ax2.set_yticks([-60,-30,0,30,60])

ax1.set_ylabel("°N")
ax2.set_xlabel("°E")

tight_layout()
savefig("maps/zfp_o3.png",dpi=200)
close(fig)
