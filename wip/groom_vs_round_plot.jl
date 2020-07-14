using Statistics
using Elefridge
using PyPlot
using StatsBase

N = 10000
m = 10000
q = 7

meanerror = fill(0f0,3,q,N)
linerror = fill(0f0,3,q,N)
logerror = fill(0f0,3,q,N)

for i in 1:N
    # sample various distributions
    for (s,D) = enumerate([ rand(Float32,m),
                            4*rand(Float32,m),
                            rand(Float32,m).+1,
                            randn(Float32,m),
                            4*randn(Float32,m),
                            randn(Float32,m).+1,
                            exp.(randn(Float32,m))])

        Dround = round(D,7)
        Dgroom = groom(D,7)
        Dshave = shave(D,7)

        # mean error
        meanerror[1,s,i] = abs(mean(D - Dround))
        meanerror[2,s,i] = abs(mean(D - Dgroom))
        meanerror[3,s,i] = abs(mean(D - Dshave))

        # max lin error
        linerror[1,s,i] = maximum(abs.(D - Dround))
        linerror[2,s,i] = maximum(abs.(D - Dgroom))
        linerror[3,s,i] = maximum(abs.(D - Dshave))

        # max log error
        logerror[1,s,i] = maximum(abs.(log10.(D./Dround)))
        logerror[2,s,i] = maximum(abs.(log10.(D./Dgroom)))
        logerror[3,s,i] = maximum(abs.(log10.(D./Dshave)))
    end
 end

## fit histogram
ion()
pygui(true)
fig,(ax1,ax2,ax3) = subplots(1,3,figsize=(10,4),sharey=true)

lwmm = 1
lw10 = 3
lw25 = 8

distr = [   "Uniform(0,1)","Uniform(0,4)","Uniform(1,2)",
            "Normal(0,1)","Normal(0,4)","Normal(1,1)",
            "Lognormal(0,1)"]

for (ax,error) in zip([ax1,ax2,ax3],[meanerror,linerror,logerror])
    for i in 1:q
        for (s,e) in enumerate([error[1,i,:],error[2,i,:],error[3,i,:]])
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
end

ax1.set_yticks([1,2,3,4,5,6,7].-0.5)
ax1.yaxis.set_tick_params(length=0)
ax2.yaxis.set_tick_params(length=0)
ax3.yaxis.set_tick_params(length=0)
ax1.set_yticklabels(distr)
ax1.set_ylim(-1.75,7)

ax1.set_xlim(1e-6,1e-2)
ax2.set_xlim(1e-3,1e0)
ax3.set_xlim(1e-3,1e0)

e = meanerror[2,1,:]
p10 = percentile(e,10)
p25 = percentile(e,25)
p50 = percentile(e,50)
p75 = percentile(e,75)
p90 = percentile(e,90)
maxx = maximum(e)
ax1.text(p10,0.1,"10%",rotation=90,ha="center",va="top",fontsize=9)
ax1.text(p25,0.1,"25%",rotation=90,ha="center",va="top",fontsize=9)
ax1.text(p50,0.1,"median",rotation=90,ha="center",va="top",fontsize=9)
ax1.text(p75,0.1,"75%",rotation=90,ha="center",va="top",fontsize=9)
ax1.text(p90,0.1,"90%",rotation=90,ha="center",va="top",fontsize=9)
ax1.text(maxx,0.1,"max",rotation=90,ha="center",va="top",fontsize=9)

ax2.plot(0,0,"C0",lw=8,label="Round-to-nearest")
ax2.plot(0,0,"C1",lw=8,label="Bit-grooming")
ax2.plot(0,0,"C2",lw=8,label="Bit-shaving")
ax2.legend(loc=4)

ax1.set_xlabel("absolute mean error")
ax2.set_xlabel("maximum absolute error")
ax3.set_xlabel("maximum decimal error")

ax1.set_title("Mean error",loc="left")
ax2.set_title("Absolute error",loc="left")
ax3.set_title("Decimal error",loc="left")

ax1.set_title("a",loc="right",fontweight="bold")
ax2.set_title("b",loc="right",fontweight="bold")
ax3.set_title("c",loc="right",fontweight="bold")


tight_layout()
