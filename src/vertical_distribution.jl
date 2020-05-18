using NetCDF, LinearAlgebra
using Statistics, StatsBase
using PyCall
pygui(:qt)
using PyPlot
using MultivariateStats
eofs = pyimport("eofs")

path = "/Users/milan/Downloads/cams"
filelist = filter(x->endswith(x,"no2.nc"),readdir(path))

ncfile = NetCDF.open(joinpath(path,filelist[1]))
no2 = ncfile.vars["no2"][:,:,:,1]
level = ncfile.vars["level"][:]
lat = ncfile.vars["latitude"][:]
lon = ncfile.vars["longitude"][:]

no2m = mean(no2,dims=(1,2))[1,1,:]
no2mi = median(no2,dims=(1,2))[1,1,:]
no2min = minimum(no2,dims=(1,2))[1,1,:]
no2max = maximum(no2,dims=(1,2))[1,1,:]
no2p10 = [percentile(vec(no2[:,:,i]),10) for i in level]
no2p90 = [percentile(vec(no2[:,:,i]),90) for i in level]

## EOF in the vertical
no2f = log10.(Matrix(reshape(no2,(:,length(level)))'))
M = fit(PCA,no2f,maxoutdim=10)
no2fhat = transform(M,no2f)
no2fr = 10f0 .^reconstruct(M,no2fhat)

##
no2rm = mean(no2fr,dims=2)[:,1]
no2rmi = median(no2fr,dims=2)[:,1]
no2rmin = minimum(no2fr,dims=2)[:,1]
no2rmax = maximum(no2fr,dims=2)[:,1]
no2rp10 = [percentile(no2fr[i,:],10) for i in level]
no2rp90 = [percentile(no2fr[i,:],90) for i in level]

##
fig,(ax,ax2) = subplots(1,2,figsize=(8,6),sharey=true,sharex=true)
ax.invert_yaxis()

ax.semilogx(no2m,level,"C0",lw=3,label="mean")
ax.plot(no2mi,level,"C0",lw=1.5,label="median")

ax.fill_betweenx(level,no2min,no2max,color="C0",alpha=0.2,label="Min-max range")
ax.fill_betweenx(level,no2p10,no2p90,color="C0",alpha=0.5,label="10-90% range")

ax2.plot(no2rm,level,"C0",lw=3,label="mean")
ax2.plot(no2rmi,level,"C0",lw=1.5,label="median")

ax2.fill_betweenx(level,no2rmin,no2rmax,color="C0",alpha=0.2,label="Min-max range")
ax2.fill_betweenx(level,no2rp10,no2rp90,color="C0",alpha=0.5,label="10-90% range")

ax.set_title(L"NO$_2$ vertical distribution")
ax2.set_title(L"NO$_2$ reconstructed from 10 PCs")
ax.set_ylabel("model level")
ax.set_xlabel("mixing ratio kg/kg")
ax2.set_xlabel("mixing ratio kg/kg")
ax.set_ylim(137,1)
ax.legend(loc=2)

tight_layout()
