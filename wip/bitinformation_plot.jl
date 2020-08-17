using PyPlot
using JLD
using Statistics, StatsBase
using ColorSchemes
using Printf

@load "/Users/milan/cams/entropy/information_all_unstructured.jld"
ICgridded = load("/Users/milan/cams/entropy/information_all.jld","IC")

# normalize
IC = mean(IC,dims=2)[:,1,:]
ICgridded = mean(ICgridded,dims=2)[:,1,:]

IC[varnames.=="t",:] .= ICgridded[varnames.=="t",:]

varnames[varnames .== "c"] .= "ch4_c"

nvars = length(varnames)

## sort into groups
aero = Array(1:15)
ozone = [35,44,45]
methane = [25,26,41]
dynamics = [33,34,54,55,53]
clouds = [22,27,28,31,32,51]
hydro = vcat(Array(36:40),46)
nitro = [42,43,52]
oopp = Array(47:50)
ces = vcat(Array(17:21),[23,24,16])
co12 = [29,30]

grouped = vcat(aero,co12,clouds,methane,ces,hydro,dynamics,nitro,ozone,oopp)
groups = [aero,co12,clouds,methane,ces,hydro,dynamics,nitro,ozone,oopp]

# sort
varnames = varnames[grouped]
IC = IC[grouped,:]

## sum over 3 dimensions
# normalise such that max entropy = 1
ICsum = copy(IC)
ICsumnan = copy(ICsum)
ICsumnan[iszero.(ICsumnan)] .= NaN


## bit start
bstart = [argmax(ICsum[i,2:end] .> 0) for i in 1:nvars]
signbits = [ICsum[i,1] .> 0 ? 1 : 0 for i in 1:nvars]

mininf = [argmin(ICsum[i,10:end])+9 for i in 1:nvars]

# manual adjustment
mininf[varnames.=="aermr08"] .= 1+8+7
mininf[varnames.=="aermr10"] .= 1+8+7
mininf[varnames.=="c5h8"] .= 1+8+4
mininf[varnames.=="ole"] .= 1+8+4

IC_noiseremoved = copy(ICsum)
for i in 1:nvars
    IC_noiseremoved[i,mininf[i]:end] .= 0
end

ICcsum = cumsum(IC_noiseremoved,dims=2)
ICcsum_norm = copy(ICcsum)
for i in 1:nvars
    ICcsum_norm[i,:] ./= ICcsum_norm[i,end]
end

inf99 = [argmax(ICcsum_norm[i,:] .> 0.99) for i in 1:nvars]

inf99x = copy(vec(hcat(inf99,inf99)'))
inf99y = copy(vec(hcat(Array(0:nvars-1),Array(1:nvars))'))

## plotting
ioff()
fig,ax1 = subplots(1,1,figsize=(8,10),sharey=true)
ax1.invert_yaxis()
tight_layout(rect=[0.05,0.08,0.93,0.98])
pos = ax1.get_position()
cax = fig.add_axes([pos.x0,0.06,pos.x1-pos.x0,0.02])

ax1right = ax1.twinx()
ax1right.invert_yaxis()

# information
cmap = ColorMap(ColorSchemes.turku.colors).reversed()
cmap_array = cmap(0:1/100000:1)
cmap_array[1,:] = [1,1,1,1]

pcm = ax1.pcolormesh(ICsumnan,vmin=0,vmax=1;cmap)
cbar = colorbar(pcm,cax=cax,orientation="horizontal")
cbar.set_label("information content [bit]")

# useful bit ranges
ax1.plot(inf99,Array(0:nvars-1),"C1",ds="steps-pre",zorder=10)

# for legend only
ax1.fill_betweenx([-1,-1],[-1,-1],[-1,-1],color="w",edgecolor="k",label="unused bits")

# grey shading
ax1.fill_betweenx(inf99y,inf99x,fill(32,length(inf99x)),alpha=0.5,color="grey",label="noise / false information")

ax1.axvline(1,color="k",lw=1,zorder=3)
ax1.axvline(9,color="k",lw=1,zorder=3)

for (ig,group) in enumerate(groups)
    y = sum([length(g) for g in groups[1:ig]])
    ax1.axhline(y,color="w",lw=1.5,zorder=2)
end

ax1.set_title("Bitwise information content",loc="left",fontweight="bold")

ax1.set_xlim(0,32)
ax1.set_ylim(nvars,0)
ax1right.set_ylim(nvars,0)

ax1.set_yticks(Array(1:nvars).-0.5)
ax1right.set_yticks(Array(1:nvars).-0.5)
ax1.set_yticklabels(varnames)
ax1right.set_yticklabels([@sprintf "%4.1f" i for i in ICcsum[:,end]])
ax1right.set_ylabel("information per value [bit]")

ax1.text(inf99[1]+0.1,0.8,"$(inf99[1]-9) significant bits",fontsize=8,color="k")
for i in 2:nvars
    ax1.text(inf99[i]+0.1,(i-1)+0.8,"$(inf99[i]-9)",fontsize=8,color="k")
end

ax1.set_xticks([1,9])
ax1.set_xticklabels([])
ax1.text(0,nvars+2,"sign",rotation=90)
ax1.text(2,nvars+1.2,"exponent")
ax1.text(10,nvars+1.2,"significant bits")

ax1.legend(loc=(0.64,0.21))

savefig("/Users/milan/git/Elefridge.jl/plots/bitinformation_all.png",dpi=100)
close(fig)
