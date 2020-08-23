using PyPlot
using JLD
using Statistics

# load linear and logarithmic quantisation
D = load("/Users/milan/cams/error/linlogroundzfp_all.jld")

Elin = D["E"][:,1,:]
Elog = D["E"][:,2,:]

@load "/Users/milan/cams/error/roundzfp_all_opt_gridded2.jld"
Eround = E[:,4,:]
Ezfp = E[:,5,:]
Ezfp[[33,34,54,55],:] = E[[33,34,54,55],6,:]

# make compression factors relative to unstructured grid
Eround[:,1] .= Eround[:,1]/(900*451)*(348528)
Ezfp[:,1] .= Ezfp[:,1]/(900*451)*(348528)


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
chydro = vcat(ces,hydro)

grouped = vcat(aero,co12,clouds,methane,chydro,dynamics,nitro,ozone,oopp)
groups = [aero,co12,clouds,methane,chydro,dynamics,nitro,ozone,oopp]

## geometric mean
cflin = 1/mean(inv,Elin[Elin[:,1].>0,1])
cflog = 1/mean(inv,Elog[Elog[:,1].>0,1])
cfround = 1/mean(inv,Eround[Eround[:,1].>0,1])
# cfround_un = 1/mean(inv,Eround_un[Eround_un[:,1].>0,1])
cfzfp = 1/mean(inv,Ezfp[Ezfp[:,1].>0,1])

## plotting
ioff()
fig,(ax1,ax2) = subplots(1,2,figsize=(10,5),sharey=true)

colours = ["grey","k","C1","C2"]
labls = ["LinQuant24","LogQuant16","Round+lossless","Zfp"]
labls2 = ["Aerosols","Carbon oxides","Clouds & water","Methane","Alkanes,\nalcohols",
            "Dynamics &\ntemperature","N&S oxides","Ozone","?"]
symblist = ["v","^",">","<","d","s","D","p","h"]

alfa = 0.7

abserri = 4     # statistic to plot, 3 = median, 4 = 90%, 5 = max
decerri = 7     # statistic to plot, 6 = median, 7 = 90%, 8 = max

# ABSOLUTE ERROR
for (i,Es) in enumerate([Elin,Elog,Eround,Ezfp])
    for (ig,g) in enumerate(groups)
        ax1.scatter(Es[g,abserri],log2.(Es[g,1]),30,colours[i],alpha=alfa,
                marker=symblist[ig],edgecolor="k")
    end
end

# DECIMAL ERROR
for (i,Es) in enumerate([Elin,Elog,Eround,Ezfp])
    for (ig,g) in enumerate(groups)
        ax2.scatter(Es[g,decerri],log2.(Es[g,1]),30,colours[i],alpha=alfa,
                    marker=symblist[ig],edgecolor="k")
    end
end

# for varname in ["d"]
#     idx = Array(1:length(varnames))[varnames .== varname]
#     for Es in [Elin,Elog,Eround,Ezfp]
#         ax1.text(Es[idx,abserri],log2.(Es[idx,1]),varname)
#         ax2.text(Es[idx,decerri],log2.(Es[idx,1]),varname)
#     end
# end

ax1.set_yticks([0,1,2,3,4,5,6,7,8,9])
ax1.set_yticks(log2.(Array(2:2:128)),minor=true)
ax1.set_yticklabels([1,2,4,8,16,32,64,128,256,512])
ax1.set_xscale("log")
ax2.set_xscale("log")

ax1.set_xlim(2e-8,9)
ax2.set_xlim(2e-7,99)

ax1.set_ylim(0,log2(128))
ax2.set_ylim(0,log2(128))

ax1.set_title("Absolute error",loc="left")
ax2.set_title("Decimal error",loc="left")

ax1.set_xlabel("norm. absolute error")
ax2.set_xlabel("decimal error")
ax1.set_ylabel("Compression factor")

ax1.set_title("a",loc="right",fontweight="bold")
ax2.set_title("b",loc="right",fontweight="bold")

for i in 1:4
    ax1.fill_betweenx([0,0],[0,0],color=colours[i],label=labls[i],alpha=0.7)
end

for i in 1:9
    ax1.scatter(0,0,30,"white",alpha=alfa,marker=symblist[i],edgecolor="k",
                    label=labls2[i])
end

ax1.legend(loc=2,scatterpoints=3,fontsize=8)

ax1.axhline(log2(cfround),color="C1",zorder=-1)
ax1.axhline(log2(cfzfp),color="C2",zorder=-1)
ax2.axhline(log2(cfround),color="C1",zorder=-1)
ax2.axhline(log2(cfzfp),color="C2",zorder=-1)

ax2.text(1e-6,log2(cfround)+0.1,"$(Int(round(cfround)))x",color="C1",fontweight="bold")
ax2.text(1e-6,log2(cfzfp)+0.1,"$(Int(round(cfzfp)))x",color="C2",fontweight="bold")

tight_layout()
savefig("/Users/milan/git/Elefridge.jl/plots/linlogroundzfp_all.png")
close(fig)
