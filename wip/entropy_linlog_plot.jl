using PyPlot
using JLD

D = load("/Users/milan/cams/entropy/gridded_linlog16.jld")

varnames = D["varnames"]
Hlin = D["Hlin"]
Hlog = D["Hlog"]

## PLOT
pygui(true)
fig,ax = subplots(1,1,figsize=(6,6))

ax.scatter(Hlin,Hlog,20,"C0",label=varnames)
#ax.legend(loc=4)

ax.plot([0,16],[0,16],"grey",lw=0.5,zorder=-1)

ax.set_xlim(0,17)
ax.set_ylim(0,17)
ax.plot([0,16],[16,16],"k")
ax.plot([16,16],[0,16],"k")
ax.text(1,16.1,"maximum entropy")

ax.set_xlabel("linear entropy [bit]")
ax.set_ylabel("logarithmic entropy [bit]")
ax.set_title("16-bit quantization entropy",loc="left")
