using PyPlot
using JLD

D = load("/Users/milan/cams/entropy/gridded_linlog16.jld")

varnames = D["varnames"]
Hlin = D["Hlin"]
Hlog = D["Hlog"]

# vargroups
aervars = [startswith(name,"aer") for name in varnames]
qvars = [name in ["clwc","crwc","cswc","ciwc"] for name in varnames]
hvars = [startswith(name,"h") for name in varnames]
humid = [name in ["q"] for name in varnames]
nitro = [startswith(name,"n") for name in varnames]
temp = [name in ["t"] for name in varnames]
lnsp = [name in ["lnsp"] for name in varnames]
co2 = [name in ["co2"] for name in varnames]
ozone = [name in ["o3","o3s","go3"] for name in varnames]
velo = [name in ["w","vo","d","z","etadot"] for name in varnames]
alc = [endswith(name,"oh") for name in varnames]
alk = [startswith(name,"c2") || endswith(name,"h8") || startswith(name,"ch4") for name in varnames ]
sul = [name in ["so2"] for name in varnames]


# rs = 10.0 .^ Array(-1:0.05:1)
# Hlintheo = [bitentropy(LinQuant16Array(exp.(r*randn(Float32,1000000)))) for r in rs]
# Hlogtheo = [bitentropy(LogQuant16Array(exp.(r*randn(Float32,1000000)))) for r in rs]

## PLOT
pygui(true)
fig,ax = subplots(1,1,figsize=(6,6))

ms = 35
al = 0.7

ax.scatter(Hlin[aervars],Hlog[aervars],ms,alpha=al,"k",label="Aerosols",marker="^")
ax.scatter(Hlin[humid],Hlog[humid],ms,alpha=al,"C4",label="Humidity",marker="o")
ax.scatter(Hlin[qvars],Hlog[qvars],ms,alpha=al,"C0",label="Cloud water species",marker="o")
ax.scatter(Hlin[hvars],Hlog[hvars],ms,alpha=al,"C2",label="Hydrogens",marker="p")
ax.scatter(Hlin[nitro],Hlog[nitro],ms,alpha=al,"C1",label="Nitrogens",marker="v")
ax.scatter(Hlin[temp],Hlog[temp],ms,alpha=al,"C3",label="Temperature",marker="<")
ax.scatter(Hlin[co2],Hlog[co2],ms,alpha=al,"C9",label="Carbon dioxide",marker=">")
ax.scatter(Hlin[lnsp],Hlog[lnsp],ms,alpha=al,"C2",label="Log surface pressure",marker="o")
ax.scatter(Hlin[ozone],Hlog[ozone],ms,alpha=al,"C5",label="Ozone",marker="d")
ax.scatter(Hlin[velo],Hlog[velo],ms,alpha=al,"C6",label="Dynamics",marker="s")
ax.scatter(Hlin[alc],Hlog[alc],ms,alpha=al,"C7",label="Alcohols",marker="h")
ax.scatter(Hlin[alk],Hlog[alk],ms,alpha=al,"C8",label="Alkanes",marker="h")
ax.scatter(Hlin[sul],Hlog[sul],ms,alpha=al,"C4",label="Sulphates",marker="h")

# ax.plot(Hlintheo,Hlogtheo)

# for (i,varname) in enumerate(varnames)
#     ax.text(Hlin[i],Hlog[i],varname,fontsize=8)
# end
ax.legend(loc=4,ncol=2,fontsize=9)

n=16

ax.plot([0,n],[0,n],"grey",lw=0.5,zorder=-1)

ax.set_xlim(0,n+1)
ax.set_ylim(0,n+1)
ax.plot([0,n],[n,n],"k")
ax.plot([n,n],[0,n],"k")
ax.text(1,n+0.1,"maximum entropy")
ax.fill_between([0,n,n,n+2],[n,n,0,0],[n+2,n+2,n+2,n+2],alpha=0.3)

ax.set_xlabel("linear packing entropy [bit]")
ax.set_ylabel("logarithmic packing entropy [bit]")
ax.set_title("$n-bit quantization entropy",loc="left")

tight_layout()
