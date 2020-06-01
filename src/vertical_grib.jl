using NetCDF, LinearAlgebra
using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
#pygui(:qt)
using PyPlot
using MultivariateStats
# eofs = pyimport("eofs")

path = "/Users/milan/Downloads/cams/ieee"
filelist = filter(x->endswith(x,"no2.grib"),readdir(path))
gribfile = xr.open_dataset(joinpath(path,filelist[1]),engine="cfgrib")
no2 = gribfile.no2.values
level = 1:size(no2)[1]

no2m = mean(no2,dims=2)[:,1]
no2mi = median(no2,dims=2)[:,1]
no2min = minimum(no2,dims=2)[:,1]
no2max = maximum(no2,dims=2)[:,1]
no2p10 = [percentile(vec(no2[i,:]),10) for i in level]
no2p90 = [percentile(vec(no2[i,:]),90) for i in level]

## EOF in the vertical
M = fit(PCA,log10.(no2),maxoutdim=50,pratio=0.9999)
no2hat = transform(M,log10.(no2))
no2r = 10f0 .^reconstruct(M,no2hat)

## reconstruction statistics
no2rm = mean(no2r,dims=2)[:,1]
no2rmi = median(no2r,dims=2)[:,1]
no2rmin = minimum(no2r,dims=2)[:,1]
no2rmax = maximum(no2r,dims=2)[:,1]
no2rp10 = [percentile(vec(no2r[i,:]),10) for i in level]
no2rp90 = [percentile(vec(no2r[i,:]),90) for i in level]

## error
decimal_error = abs.(log10.(no2 ./ no2r))

dem = mean(decimal_error,dims=2)[:,1]
demi = median(decimal_error,dims=2)[:,1]
demin = minimum(decimal_error,dims=2)[:,1]
demax = maximum(decimal_error,dims=2)[:,1]
dep10 = [percentile(vec(decimal_error[i,:]),10) for i in level]
dep90 = [percentile(vec(decimal_error[i,:]),90) for i in level]

## highest loadings
P = projection(M)
Ps = zeros(Int64,size(P))
for i in 1:size(Ps)[2]
    Ps[:,i] = sortperm(abs.(P[:,i]))
end
L = zeros(Int64,size(P)[2])
L[1] = Ps[end,1]
for i in 2:length(L)
    j = 0
    while Ps[end-j,i] in L
        j += 1
    end
    L[i] = Ps[end-j,i]
end

## plot
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

for i in level
    xvals = no2[i,(decimal_error[i,:] .>= 0.5f0) .* (decimal_error[i,:] .< 1f0)]
    ax.scatter(xvals,i*ones(length(xvals)),10,"yellow",alpha=0.05)

    xvals = no2[i,(decimal_error[i,:] .>= 1f0) .* (decimal_error[i,:] .< 2f0)]
    ax.scatter(xvals,i*ones(length(xvals)),10,"C1",alpha=0.1)
end

ax.scatter([1,2,3],[0,0,0],10,"yellow",alpha=0.5,label="error > 0.5")
ax.scatter([1,2,3],[0,0,0],10,"C1",alpha=0.5,label="error > 1")

ax3.semilogx(dem,level,"C0",lw=3,label="mean")
ax3.plot(demi,level,"C0",lw=1.5,label="median")
ax3.fill_betweenx(level,demin,demax,color="C0",alpha=0.2,label="Min-max range")
ax3.fill_betweenx(level,dep10,dep90,color="C0",alpha=0.5,label="10-90% range")

cmap = ColorMap("viridis")
for i in 1:length(L)
    ax2.plot([1e-14,1e-13],[L[i],L[i]].-0.5,color=cmap(i/length(L)),lw=2)
    ax2.text(2e-13,L[i],"$i")
end
ax2.text(1.5e-14,8,"modes")

ax.set_title(L"NO$_2$ vertical distribution")
ax2.set_title(L"NO$_2$ reconstructed from 50 PCs")
ax3.set_title("L1 log-error")
ax.set_ylabel("model level")
ax.set_xlabel("mixing ratio kg/kg")
ax2.set_xlabel("mixing ratio kg/kg")
ax.set_ylim(137,1)
ax.legend(loc=2,scatterpoints=3)

ax.set_xlim(1e-14,1e-6)
ax2.set_xlim(1e-14,1e-6)

tight_layout()
