using PyPlot
using JLD
<<<<<<< HEAD
using Statistics
=======
>>>>>>> 5baed54df0982b9ab5d73d7af6a9f125ceb1e66c

@load "/Users/milan/cams/error/linlogroundzfp_all.jld"

sort_out = .~[varname in ["w","etadot","d","vo"] for varname in varnames ]
varnames = varnames[sort_out]
E = E[sort_out,:,:]

## compression level
complev = fill(3,length(varnames))
complev[varnames .== "t"] .= 2
complev[varnames .== "co2"] .= 1
complev[varnames .== "c"] .= 2
complev[varnames .== "aermr18"] .= 2

complevz = fill(3,length(varnames))
complevz[varnames .== "t"] .= 2
complevz[varnames .== "co2"] .= 2
complevz[varnames .== "aermr18"] .= 2


Eround = fill(0f0,length(varnames),length(stats))
Ezfp = fill(0f0,length(varnames),length(stats))

for i in 1:length(varnames)
    Eround[i,:] = E[i,5+complev[i],:]
    Ezfp[i,:] = E[i,8+complevz[i],:]
end

Elin = E[:,1,:]
Elog = E[:,2,:]

## geometric mean
cflin = 1/mean(inv,Elin[:,1])
cflog = 1/mean(inv,Elog[:,1])
cfround = 1/mean(inv,Eround[:,1])
cfzfp = 1/mean(inv,Ezfp[:,1])

## plotting
ioff()
fig,(ax1,ax2) = subplots(1,2,figsize=(10,5),sharey=true)

colours = ["grey","k","C1","C2"]
labls = ["LinQuant24","LogQuant16","Round+lossless","Zfp"]
alfa = 0.7

# ABSOLUTE ERROR
for (i,Es) in enumerate([Elin,Elog,Eround,Ezfp])
    ax1.scatter(Es[:,4],log2.(Es[:,1]),30,colours[i],alpha=alfa,
                edgecolor="k",label=labls[i])
end

# DECIMAL ERROR
for (i,Es) in enumerate([Elin,Elog,Eround,Ezfp])
    ax2.scatter(Es[:,7],log2.(Es[:,1]),30,colours[i],alpha=alfa,
                edgecolor="k")
end

# for (i,varname) in enumerate(varnames)
#     ax2.text(Ezfp[i,7],log2.(Ezfp[i,1]),varname)
# end
#
# for (i,varname) in enumerate(varnames)
#     ax1.text(Eround[i,7],log2.(Eround[i,1]),varname)
# end

ax1.set_yticks([0,1,2,3,4,5,6,7,8,9])
ax1.set_yticks(log2.(Array(2:2:100)),minor=true)
ax1.set_yticklabels([1,2,4,8,16,32,64,128,256,512])
ax1.set_xscale("log")
ax2.set_xscale("log")

ax1.set_xlim(1e-8,0.9)
ax2.set_xlim(2e-7,20)

ax1.set_ylim(0,log2(100))
ax2.set_ylim(0,log2(100))

ax1.set_title("Absolute error",loc="left")
ax2.set_title("Decimal error",loc="left")

ax1.set_xlabel("norm. absolute error")
ax2.set_xlabel("decimal error")
ax1.set_ylabel("Compression factor")

ax1.set_title("a",loc="right",fontweight="bold")
ax2.set_title("b",loc="right",fontweight="bold")

ax1.legend(loc=2,scatterpoints=3)

ax1.axhline(log2(cfround),color="C1",zorder=-1)
ax1.axhline(log2(cfzfp),color="C2",zorder=-1)
ax2.axhline(log2(cfround),color="C1",zorder=-1)
ax2.axhline(log2(cfzfp),color="C2",zorder=-1)

ax2.text(1e-6,log2(cfround)+0.1,"$(Int(round(cfround)))x",color="C1",fontweight="bold")
ax2.text(1e-6,log2(cfzfp)+0.1,"$(Int(round(cfzfp)))x",color="C2",fontweight="bold")

tight_layout()
savefig("/Users/milan/git/Elefridge.jl/plots/linlogroundzfp_all.png")
close(fig)
