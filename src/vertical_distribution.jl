using NetCDF, LinearAlgebra
using Statistics, StatsBase
using PyCall
pygui(:qt)
using PyPlot
using MultivariateStats
# eofs = pyimport("eofs")

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
M = fit(PCA,no2f,maxoutdim=20,pratio=0.999)
no2fhat = transform(M,no2f)
no2fr = 10f0 .^reconstruct(M,no2fhat)

## reconstruction statistics
no2rm = mean(no2fr,dims=2)[:,1]
no2rmi = median(no2fr,dims=2)[:,1]
no2rmin = minimum(no2fr,dims=2)[:,1]
no2rmax = maximum(no2fr,dims=2)[:,1]
no2rp10 = [percentile(no2fr[i,:],10) for i in level]
no2rp90 = [percentile(no2fr[i,:],90) for i in level]

## error
decimal_error = abs.(log10.(reshape(no2,(:,length(level)))' ./ no2fr))

dem = mean(decimal_error,dims=2)[:,1]
demi = median(decimal_error,dims=2)[:,1]
demin = minimum(decimal_error,dims=2)[:,1]
demax = maximum(decimal_error,dims=2)[:,1]
dep10 = [percentile(decimal_error[i,:],10) for i in level]
dep90 = [percentile(decimal_error[i,:],90) for i in level]

## by number of modes
M = fit(PCA,no2f,maxoutdim=80,pratio=0.999999)
no2fhat = transform(M,no2f)
P = projection(M)
n = 71
gm = Array{Float32,2}(undef,n,6)

# for i in 1:size(no2f)[2]
    # no2f[:,i] -= mean(M)
# end

for i in 1:n
    println(i)
    no2fr = P[:,1:9+i]*no2fhat[1:9+i,:]
    decimal_error = abs.(no2f-no2fr)
    gm[i,1] = mean(decimal_error)
    gm[i,2] = median(decimal_error)
    gm[i,3] = minimum(decimal_error)
    gm[i,4] = maximum(decimal_error)
    gm[i,5] = percentile(vec(decimal_error),10)
    gm[i,6] = percentile(vec(decimal_error),90)
end

##
fig,(ax,ax2,ax3) = subplots(1,3,figsize=(10,6),sharey=true)
ax.invert_yaxis()

ax.semilogx(no2m,level,"C0",lw=3,label="mean")
ax.plot(no2mi,level,"C0",lw=1.5,label="median")
ax.fill_betweenx(level,no2min,no2max,color="C0",alpha=0.2,label="Min-max range")
ax.fill_betweenx(level,no2p10,no2p90,color="C0",alpha=0.5,label="10-90% range")

ax2.semilogx(no2rm,level,"C0",lw=3,label="mean")
ax2.plot(no2rmi,level,"C0",lw=1.5,label="median")
ax2.fill_betweenx(level,no2rmin,no2rmax,color="C0",alpha=0.2,label="Min-max range")
ax2.fill_betweenx(level,no2rp10,no2rp90,color="C0",alpha=0.5,label="10-90% range")

ax3.semilogx(dem,level,"C0",lw=3,label="mean")
ax3.plot(demi,level,"C0",lw=1.5,label="median")
ax3.fill_betweenx(level,demin,demax,color="C0",alpha=0.2,label="Min-max range")
ax3.fill_betweenx(level,dep10,dep90,color="C0",alpha=0.5,label="10-90% range")

ax.set_title(L"NO$_2$ vertical distribution")
ax2.set_title(L"NO$_2$ reconstructed from 20 PCs")
ax3.set_title("L1 log-error")
ax.set_ylabel("model level")
ax.set_xlabel("mixing ratio kg/kg")
ax2.set_xlabel("mixing ratio kg/kg")
ax.set_ylim(137,1)
ax.legend(loc=2)

ax.set_xlim(1e-15,1e-7)
ax2.set_xlim(1e-15,1e-7)

tight_layout()

##
modes = Array(10:79)

fig,ax1 = subplots(1,1)

ax1.semilogy(modes,gm[:,1],"C0",lw=3,label="mean")
ax1.plot(modes,gm[:,2],"C0",lw=1.5,label="median")
ax1.fill_between(modes,gm[:,3],gm[:,4],color="C0",alpha=0.2,label="Min-max range")
ax1.fill_between(modes,gm[:,5],gm[:,6],color="C0",alpha=0.5,label="10-90% range")
ax1.legend(loc=1)

ax1.set_title("Decimal error by modes",loc="left")
ax1.set_ylabel("Decimal error")
ax1.set_xlabel("number of modes used")
ax1.set_xlim(10,79)
tight_layout()
