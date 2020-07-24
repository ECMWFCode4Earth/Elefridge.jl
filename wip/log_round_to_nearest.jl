#using Elefridge
using PyPlot

nbits = 8
nmax = 7
mi,ma = 1e-20,1e20
logmin = log(mi)
logmax = log(ma)
dlog = (logmax-logmin)/(2^nbits-2)              # 2^n-1-1 as 0x00 is reserved for 0
quants = exp.(logmin .+ Array(0:nmax)*dlog)     # the first quantised numbers
half_quants = quants[1:end-1] + diff(quants)/2         # half-way in lin-space
halflog_quants = exp.(logmin .+ (Array(0:nmax-1).+.5)*dlog)    # half-way in log-space
dlin = quants[2]-quants[1]

A = Array(quants[1]:dlin/100:quants[end])

# round to nearest in log space
L1 = round.((log.(A).-logmin)/dlog).+1
A1 = exp.(logmin .+ (L1.-1)*dlog)

B = vcat(A,ma)
A0 = Array(LogQuant8Array(B))

# round to nearest in lin space
c = 1/2 - 1/dlog*log(mi*(exp(dlog)+1)/2)

L2 = round.(c.+log.(A)/dlog).+1
A2 = exp.(logmin .+ (L2.-1)*dlog)


## plot
pygui(true)
fig,ax = subplots(1,1)
ax.plot(A,A1,"k--",lw=2,label="Round-to-nearest in log-space")
ax.plot(A,A2,"C1--",lw=2,label="Round-to-nearest in lin-space")
ax.plot(B,A0,"C2.-")

ax.plot([quants[1],half_quants[1]],[quants[1],quants[1]],"C1",lw=3.5)
ax.plot([quants[1],halflog_quants[1]],[quants[1],quants[1]],"k",lw=1.5)
for i in 1:length(halflog_quants)-1
    x0 = half_quants[i]
    x1 = half_quants[i+1]
    y = quants[i+1]
    ax.plot([x0,x1],[y,y],"C1",lw=3.5)
    x0 = halflog_quants[i]
    x1 = halflog_quants[i+1]
    ax.plot([x0,x1],[y,y],"k",lw=1.5)
end
ax.plot([half_quants[end],quants[end]],[quants[end],quants[end]],"C1",lw=3.5)
ax.plot([halflog_quants[end],quants[end]],[quants[end],quants[end]],"k",lw=1.5)

for (i,x) in enumerate(halflog_quants)
    labl = i==1 ? "Halfway in log-space" : nothing
    ax.axvline(x,color="k",lw=1,label=labl)
end

for (i,x) in enumerate(half_quants)
    labl = i==1 ? "Halfway in lin-space" : nothing
    ax.axvline(x,color="C1",lw=1,label=labl)
end

ax.set_yticks(quants)
ax.set_yticklabels([latexstring("x_{$i}") for i in 1:length(quants)])
ax.set_xticks(quants)
ax.set_xticklabels([latexstring("x_{$i}") for i in 1:length(quants)])
ax.set_xlabel(L"$x$ at full precision")
ax.set_ylabel(L"quantized $x$")

ax.set_xlim(quants[1],quants[end])
ax.set_ylim(quants[1]-dlin/2,quants[end]+dlin)

ax.legend(loc=2,framealpha=1)
ax.set_title("Round-to-nearest for log quantization",loc="left")

tight_layout()
