using JLD
using PyPlot
@load "/Users/milan/cams/error/linvslog_unstructured.jld"

nmethods,nvars,nstats = size(linerror)

# sort by groups
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

allvars = allvars[reverse(grouped)]
meanerror = meanerror[:,reverse(grouped)]
linerror = linerror[:,reverse(grouped),:]
logerror = logerror[:,reverse(grouped),:]

allvars[allvars .== "c"] .= "ch4_c"

## fit histogram
ioff()
fig,(ax1,ax2,ax3) = subplots(1,3,figsize=(10,8),sharey=true)

lwmm = 1
lw10 = 2
lw25 = 3

colours = ["C0","k","C1","C3"]

for i in 1:nvars

    # MEAN ERROR
    for j in 1:nmethods
        me = meanerror[j,i]
        y = i-j/5
        ax1.set_xscale("log")
        ax1.scatter(me,y,30,colours[j],edgecolor="k",zorder=5)
        ax1.axhline(i-1,color="grey",lw=0.5)
    end

    # ABSOLUTE ERROR
    for j in 1:nmethods
        mi,p10,p25,p50,p75,p90,ma = linerror[j,i,:]
        y = i-j/5
        ax2.semilogx([p50,ma],[y,y],colours[j],lw=lwmm)
        ax2.semilogx([p50,p90],[y,y],colours[j],lw=lw10)
        ax2.scatter(p50,y,30,colours[j],edgecolor="k",zorder=5)
        ax2.axhline(i-1,color="grey",lw=0.5)
    end

    # DECIMAL ERROR
    for j in 1:nmethods
        mi,p10,p25,p50,p75,p90,ma = logerror[j,i,:]
        y = i-j/5
        ax3.semilogx([p50,ma],[y,y],colours[j],lw=lwmm)
        ax3.semilogx([p50,p90],[y,y],colours[j],lw=lw10)
        ax3.scatter(p50,y,30,colours[j],edgecolor="k",zorder=5)
        ax3.axhline(i-1,color="grey",lw=0.5)
    end

end

for ax in (ax1,ax2,ax3)
    for x in [4,7,10,15,32,38,40]
        ax.axhline(x,color="grey",lw=2)
    end
end

ax1.set_yticks(Array(1:nvars).-0.5)
ax1.yaxis.set_tick_params(length=0)
ax2.yaxis.set_tick_params(length=0)
ax3.yaxis.set_tick_params(length=0)
ax1.set_yticklabels(allvars)
ax1.set_ylim(-6.5,nvars+0.25)

ax1.set_xlim(1e-12,9e-1)
ax2.set_xlim(2e-22,1e1)
ax3.set_xlim(1.1e-8,1e2)

alfa=0.25
for i in 1:nvars
    if meanerror[1,i] >= meanerror[3,i] && meanerror[2,i] < meanerror[4,i]
        ax1.fill_between(ax1.get_xlim(),[i-1,i-1],[i,i],color="C9",alpha=alfa)
    end

    for (ax,e) in zip([ax2,ax3],[linerror,logerror])
        if e[1,i,4] >= e[3,i,4] && e[2,i,4] < e[4,i,4]
            ax.fill_between(ax.get_xlim(),[i-1,i-1],[i,i],color="C9",alpha=alfa)
        end
    end
end

for i in 1:nvars
    if meanerror[2,i] >= meanerror[4,i]
        ax1.fill_between(ax1.get_xlim(),[i-1,i-1],[i,i],color="C2",alpha=alfa)
    end

    for (ax,e) in zip([ax2,ax3],[linerror,logerror])
        if e[2,i,4] >= e[4,i,4]
            ax.fill_between(ax.get_xlim(),[i-1,i-1],[i,i],color="C2",alpha=alfa)
        end
    end
end

mi,p10,p25,p50,p75,p90,ma = linerror[4,1,:]
ax2.text(p50,-0.5,"median",rotation=90,ha="center",va="top",fontsize=9)
ax2.text(p90,-0.5,"90%",rotation=90,ha="center",va="top",fontsize=9)
ax2.text(ma,-0.5,"max",rotation=90,ha="center",va="top",fontsize=9)

ax1.scatter(0,0,30,"C0",edgecolor="k",label="LinQuant16")
ax1.scatter(0,0,30,"k",edgecolor="k",label="LinQuant24")
ax1.scatter(0,0,30,"C1",edgecolor="k",label="LogQuant16")
ax1.scatter(0,0,30,"C3",edgecolor="k",label="LogQuant24")
ax3.fill_between(ax3.get_xlim(),[-10,-10],[-9,-9],color="C9",alpha=alfa,label="at 16 bit")
ax3.fill_between(ax3.get_xlim(),[-10,-10],[-9,-9],color="C2",alpha=alfa,label="at both 16 and 24 bit")
ax1.legend(loc=4,ncol=2,fontsize=8,scatterpoints=3)
ax3.legend(loc=4,fontsize=9,title="Log better than LinQuantization")

ax1.set_xlabel("norm. mean error")
ax2.set_xlabel("norm. absolute error")
ax3.set_xlabel("decimal error")

ax1.set_title("Mean error",loc="left")
ax2.set_title("Absolute error",loc="left")
ax3.set_title("Decimal error",loc="left")

ax1.set_title("a",loc="right",fontweight="bold")
ax2.set_title("b",loc="right",fontweight="bold")
ax3.set_title("c",loc="right",fontweight="bold")

tight_layout()
savefig("/Users/milan/git/Elefridge.jl/plots/linvslog_sorted.png",dpi=200)
close(fig)
