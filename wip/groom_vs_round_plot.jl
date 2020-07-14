using Statistics
using Elefridge
using PyPlot
using StatsBase

N = 10000
m = 10000
q = 7

errorsr = fill(0f0,q,N)
errorsg = fill(0f0,q,N)
errorss = fill(0f0,q,N)

for i in 1:N
    # sample various distributions
    for (s,D) = enumerate([ rand(Float32,m),
                            4*rand(Float32,m),
                            rand(Float32,m).+1,
                            randn(Float32,m),
                            4*randn(Float32,m),
                            randn(Float32,m).+1,
                            exp.(randn(Float32,m))])
        errorsr[s,i] = mean(D - round(D,7))
        errorsg[s,i] = mean(D - groom(D,7))
        errorss[s,i] = mean(D - shave(D,7))
    end
 end

## fit histogram
ion()
pygui(true)
fig,ax = subplots(1,1)

lwmm = 1
lw10 = 3
lw25 = 8

distr = [   "Uniform(0,1)","Uniform(0,4)","Uniform(1,2)",
            "Normal(0,1)","Normal(0,4)","Normal(1,1)",
            "Lognormal(0,1)"]

for i in 1:q
    eround = abs.(errorsr[i,:])
    egroom = abs.(errorsg[i,:])
    eshave = abs.(errorss[i,:])

    for (s,e) in enumerate([eround,egroom,eshave])
        y = i-s/4
        c = "C"*string(s-1)
        ax.semilogx([minimum(e),maximum(e)],[y,y],c,lw=lwmm)
        p10 = percentile(e,10)
        p25 = percentile(e,25)
        p50 = percentile(e,50)
        p75 = percentile(e,75)
        p90 = percentile(e,90)
        ax.semilogx([p10,p90],[y,y],c,lw=lw10)
        ax.semilogx([p25,p75],[y,y],c,lw=lw25)
        ax.scatter(p50,y,100,c,edgecolor="k",zorder=5)
        ax.axhline(i,color="k",lw=0.5)
    end
end

ax.set_yticks([1,2,3,4,5,6,7].-0.5)
ax.yaxis.set_tick_params(length=0)
ax.set_yticklabels(distr)

ax.set_xlim(1e-6,1e-2)
ax.set_ylim(-1.75,7)

e = abs.(errorsg[1,:])
p10 = percentile(e,10)
p25 = percentile(e,25)
p50 = percentile(e,50)
p75 = percentile(e,75)
p90 = percentile(e,90)
maxx = maximum(e)
ax.text(p10,0.1,"10%",rotation=90,ha="center",va="top",fontsize=10)
ax.text(p25,0.1,"25%",rotation=90,ha="center",va="top",fontsize=10)
ax.text(p50,0.1,"median",rotation=90,ha="center",va="top",fontsize=10)
ax.text(p75,0.1,"75%",rotation=90,ha="center",va="top",fontsize=10)
ax.text(maxx,0.1,"max",rotation=90,ha="center",va="top",fontsize=10)

ax.plot(0,0,"C0",lw=8,label="Round-to-nearest")
ax.plot(0,0,"C1",lw=8,label="Bit-grooming")
ax.plot(0,0,"C2",lw=8,label="Bit-shaving")
ax.legend(loc=4)

# ax.set_yticks([])
ax.set_xlabel("absolute mean error")
# ax.set_ylabel("N")
ax.set_title("Bias comparison",loc="left")
tight_layout()
