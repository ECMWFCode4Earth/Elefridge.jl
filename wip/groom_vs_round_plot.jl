using Statistics
using Elefridge
using PyPlot
using StatsBase

N = 10000
m = 100000
q = 7

errorsr = fill(0f0,q,N)
errorsg = fill(0f0,q,N)

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
    end
 end

## fit histogram
ion()
pygui(true)
fig,ax = subplots(1,1)

lwmm = 1
lw10 = 3
lw25 = 8

distr = [   "Uniform(0,1)","Uniform(1,2)","Uniform(0,4)",
            "Normal(0,1)","Normal(1,1)","Normal(0,4)",
            "Lognormal(0,1)"]

for i in 1:q
    eround = abs.(errorsr[i,:])
    egroom = abs.(errorsg[i,:])

    for (s,e) in enumerate([eround,egroom])
        y = i-s/3
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
    end
    ax.text(1.5e-4,i-1/2,distr[i],va="center")
end


ax.set_xlim(1e-7,1e-3)
ax.set_ylim(-1,7)

e = abs.(errorsg[1,:])
p10 = percentile(e,10)
p25 = percentile(e,25)
p50 = percentile(e,50)
p75 = percentile(e,75)
p90 = percentile(e,90)
maxx = maximum(e)
ax.text(p10,0.1,"10%",rotation=90,ha="center",va="top",fontsize=9)
ax.text(p25,0.1,"25%",rotation=90,ha="center",va="top",fontsize=9)
ax.text(p50,0.1,"median",rotation=90,ha="center",va="top",fontsize=9)
ax.text(p75,0.1,"75%",rotation=90,ha="center",va="top",fontsize=9)
ax.text(maxx,0.1,"max",rotation=90,ha="center",va="top",fontsize=9)

ax.plot(0,0,"C0",lw=8,label="Round-to-nearest")
ax.plot(0,0,"C1",lw=8,label="Bit-grooming")
ax.legend(loc=4)

ax.set_yticks([])
ax.set_xlabel("absolute mean error")
# ax.set_ylabel("N")
ax.set_title("Bias: bit-grooming vs round-to-nearest",loc="left")
tight_layout()
