using PyPlot
using JLD
using Statistics, StatsBase
using ColorSchemes
using Printf

@load "/Users/milan/cams/entropy/information_all.jld"

nvars = length(varnames)

# sum over 3 dimensions
# normalise such that max entropy = 1
ICsum = sum(IC,dims=2)[:,1,:]/3
ICsumsum = sum(IC,dims=(2,3))[:,1,1]/3

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


## plotting
ioff()
fig,ax1 = subplots(1,1,figsize=(6,8),sharey=true)
ax1.invert_yaxis()
tight_layout(rect=[0.08,0.1,0.95,0.98])
pos = ax1.get_position()
cax = fig.add_axes([pos.x0,0.06,pos.x1-pos.x0,0.02])

ax1right = ax1.twinx()
ax1right.invert_yaxis()

# information
cmap = ColorMap(ColorSchemes.turku.colors).reversed()
pcm = ax1.pcolormesh(ICsumnan,vmin=0,vmax=1;cmap)
cbar = colorbar(pcm,cax=cax,orientation="horizontal")
cbar.set_label("[bit]")

# useful bit ranges
# ax1.plot(signbits,Array(0:nvars-1),"C1",ds="steps-pre",zorder=10)
# ax1.plot(bstart,Array(0:nvars-1),"C1",ds="steps-pre",zorder=10)
ax1.plot(inf99,Array(0:nvars-1),"C1",ds="steps-pre",zorder=10)
# ax1.plot(mininf,Array(0:nvars-1),"C2",ds="steps-pre",zorder=10)

ax1.fill_betweenx(Array(0:nvars-1),inf99,fill(32,nvars),alpha=0.5,color="grey")

ax1.axvline(1,color="k",lw=1)
# ax1.axvline(1,color="grey",lw=2,ls="--")
ax1.axvline(9,color="k",lw=1)
# ax1.axvline(9,color="grey",lw=2,ls="--")

ax1.set_title("Bitwise information content",loc="left",fontweight="bold")

ax1.set_xlim(0,32)
ax1.set_ylim(nvars,0)
ax1right.set_ylim(nvars,0)

ax1.set_yticks(Array(1:nvars).-0.5)
ax1right.set_yticks(Array(1:nvars).-0.5)
ax1.set_yticklabels(varnames)
ax1right.set_yticklabels([@sprintf "%4.1fbit" i for i in ICsumsum])

ax1.set_xticks([1,9])
ax1.set_xticklabels([])
ax1.text(0,nvars+3.5,"sign",rotation=90)
ax1.text(2,nvars+2,"exponent")
ax1.text(10,nvars+2,"significant bits")

savefig("/Users/milan/git/Elefridge.jl/plots/bitinformation_all.png")
close(fig)
