# using Statistics, StatsBase
# using PyCall
# xr = pyimport("xarray")
# using PyPlot
#
# path = "/Users/milan/cams"
# allvars = ["no2","go3","so2","aermr04","aermr05","aermr06","ch4","co"]
#
# M = Array{Float32,3}(undef,length(allvars),6,137)
#
# for (i,vari) in enumerate(allvars)
#     filelist = filter(x->endswith(x,vari*".grib"),readdir(path))
#     X = xr.open_dataset(joinpath(path,filelist[1]),engine="cfgrib")[vari].data
#
#     M[i,1,:] = mean(X,dims=2)[:,1]
#     M[i,2,:] = median(X,dims=2)[:,1]
#     M[i,3,:] = minimum(X,dims=2)[:,1]
#     M[i,4,:] = maximum(X,dims=2)[:,1]
#     M[i,5,:] = [percentile(vec(X[i,:]),10) for i in level]
#     M[i,6,:] = [percentile(vec(X[i,:]),90) for i in level]
# end

# plotting
level = 1:137
alphabet = ["a","b","c","d","e","f","g","h","i","j"]

fig,axs = subplots(2,4,figsize=(10,8),sharey=true)
axs[1,1].invert_yaxis()

for (i,ax) in enumerate(axs)
    ax.semilogx(M[i,1,:],level,"C0",lw=3,label="mean")
    ax.plot(M[i,2,:],level,"C0",lw=1.5,label="median")
    ax.fill_betweenx(level,M[i,3,:],M[i,4,:],color="C0",alpha=0.2,label="Min-max range")
    ax.fill_betweenx(level,M[i,5,:],M[i,6,:],color="C0",alpha=0.5,label="10-90% range")
    ax.set_title(allvars[i],loc="left")
    ax.set_title(alphabet[i],loc="right",fontweight="bold")
end

axs[2,1].legend(loc=2)

axs[1,1].set_ylim(137,1)
axs[1,1].set_ylabel("model level")
axs[2,1].set_ylabel("model level")
axs[2,1].set_xlabel("values")
axs[2,2].set_xlabel("values")
axs[2,3].set_xlabel("values")
axs[2,4].set_xlabel("values")

tight_layout()
