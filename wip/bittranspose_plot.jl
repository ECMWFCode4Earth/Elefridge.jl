using PyPlot
using JLD
using StatsBase, Statistics

@load "/Users/milan/cams/entropy/compression_transpose_all_zstd_1_8.jld"

ZstdL3 = C[:,1,:] ./ C[:,2,:]
ZstdL10 = C[:,3,:] ./ C[:,4,:]
ZstdL22 = C[:,5,:] ./ C[:,6,:]

nready = argmin(C[:,1,1])-1

# geometric mean
ZstdL3m = (mean(inv,C[1:nready,2,:],dims=1) ./ mean(inv,C[1:nready,1,:],dims=1))[1,:]
ZstdL10m = (mean(inv,C[1:nready,4,:],dims=1) ./ mean(inv,C[1:nready,3,:],dims=1))[1,:]
ZstdL22m = (mean(inv,C[1:nready,6,:],dims=1) ./ mean(inv,C[1:nready,5,:],dims=1))[1,:]

# geometric mean
cfL3m = 1 ./ mean(inv,C[1:nready,2,:],dims=1)[1,:]
cfL10m = 1 ./ mean(inv,C[1:nready,4,:],dims=1)[1,:]
cfL22m = 1 ./ mean(inv,C[1:nready,6,:],dims=1)[1,:]

# percentiles
ZstdL3a = [percentile(C[1:nready,1,s]./C[1:nready,2,s],10) for s in 1:length(sbits)]
ZstdL3b = [percentile(C[1:nready,1,s]./C[1:nready,2,s],90) for s in 1:length(sbits)]

ZstdL10a = [percentile(C[1:nready,3,s]./C[1:nready,4,s],10) for s in 1:length(sbits)]
ZstdL10b = [percentile(C[1:nready,3,s]./C[1:nready,4,s],90) for s in 1:length(sbits)]

ZstdL22a = [percentile(C[1:nready,5,s]./C[1:nready,6,s],10) for s in 1:length(sbits)]
ZstdL22b = [percentile(C[1:nready,5,s]./C[1:nready,6,s],90) for s in 1:length(sbits)]

cfL3a = [percentile(C[1:nready,1,s],10) for s in 1:length(sbits)]
cfL3b = [percentile(C[1:nready,1,s],90) for s in 1:length(sbits)]

cfL10a = [percentile(C[1:nready,3,s],10) for s in 1:length(sbits)]
cfL10b = [percentile(C[1:nready,3,s],90) for s in 1:length(sbits)]

cfL22a = [percentile(C[1:nready,5,s],10) for s in 1:length(sbits)]
cfL22b = [percentile(C[1:nready,5,s],90) for s in 1:length(sbits)]


## PLOT
fig,(ax2,ax) = subplots(2,1,figsize=(8,8),sharex=true)

colours = ["darkslategray","forestgreen","goldenrod"]

ax.fill_between(sbits,ZstdL22a,ZstdL22b,color=colours[3],alpha=0.25)
ax.fill_between(sbits,ZstdL10a,ZstdL10b,color=colours[2],alpha=0.25)
ax.fill_between(sbits,ZstdL3a,ZstdL3b,color=colours[1],alpha=0.25)

ax.plot(sbits,ZstdL22m,colours[3],lw=2.5,label="High compression (ZstdL22)")
ax.plot(sbits,ZstdL10m,colours[2],lw=2.5,label="Medium compression (ZstdL10)")
ax.plot(sbits,ZstdL3m,colours[1],lw=2.5,label="Low compression (ZstdL3)")

ax.text(19,1.5,"No transpose\nbetter",rotation=0)
ax.text(19,0.2,"Bit transpose\nbetter",rotation=0)

ax.axhline(1,color="k",ls="--")


ax2.fill_between(sbits,log2.(cfL22a),log2.(cfL22b),color=colours[3],alpha=0.25)
ax2.fill_between(sbits,log2.(cfL10a),log2.(cfL10b),color=colours[2],alpha=0.25)
ax2.fill_between(sbits,log2.(cfL3a),log2.(cfL3b),color=colours[1],alpha=0.25)

ax2.plot(sbits,log2.(cfL22m),colours[3],lw=2.5,label="High compression (ZstdL22)")
ax2.plot(sbits,log2.(cfL10m),colours[2],lw=2.5,label="Medium compression (ZstdL10)")
ax2.plot(sbits,log2.(cfL3m),colours[1],lw=2.5,label="Low compression (ZstdL3)")

ax.set_xticks(sbits)

ax.set_xlim(sbits[1],sbits[end])
ax.set_ylim(0,3)
ax2.set_ylim(0,6)

ax2.set_yticks([0,1,2,3,4,5,6])
ax2.set_yticks(log2.(Array(2:2:128)),minor=true)
ax2.set_yticklabels([1,2,4,8,16,32,64,128])

ax.set_title("Bittranspose's impact on compression factors",loc="left")
ax2.set_title("Compression factors",loc="left")

ax.set_title("b",loc="right",fontweight="bold")
ax2.set_title("a",loc="right",fontweight="bold")

ax2.legend(loc=1)

ax2.set_ylabel("Compression factor")
ax.set_ylabel("Compression factor ratio")
ax.set_xlabel("Significant bits retained")

tight_layout()
savefig("/Users/milan/git/Elefridge.jl/plots/bittranpose_zstd.png",dpi=100)
close(fig)
