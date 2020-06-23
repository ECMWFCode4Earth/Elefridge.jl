using PyPlot
using JLD

path = "/Users/milan/cams"
allvars = ["no2","go3","so2","aermr04","aermr05","aermr06","co","no"]
allvarsLN = [L"NO$_2$",L"O$_3$",L"SO$_2$","AERMR04","AERMR05","AERMR06","CO","NO"]

D = load("cams/error/$vari.jld")
EV = D["EV"]
lz4hcc = D["lz4"]

EVs = Array{Float32,3}(undef,length(allvars),size(EV)[1],size(EV)[2])
LZs = Array{Float32,1}(undef,length(allvars))

for (i,vari) in enumerate(allvars)
    D = load("cams/error/$vari.jld")
    EVs[i,:,:] = D["EV"]
    LZs[i] = D["lz4"]
end

## theoretical
theosize = [32/i for i in 32:-1:1]
x0 = 1f0
theoerror = [log10(nextfloat(x0,2^i)) for i in 0:31]

## plotting
pygui(true)

cm1 = ColorMap("viridis")
cm2 = ColorMap("inferno")

function RGB(cv::Tuple)
    r = UInt8(round(cv[1]*255))
    g = UInt8(round(cv[2]*255))
    b = UInt8(round(cv[3]*255))
    return "#"*repr(r)[3:4]*repr(g)[3:4]*repr(b)[3:4]
end

colours = ["0.6","0.25",RGB(cm2(0.8)),RGB(cm2(0.6)),RGB(cm2(0.4)),RGB(cm2(0.9)),
                    RGB(cm1(0.6)),RGB(cm1(0.4)),RGB(cm1(0.2)),RGB(cm1(0.8))]

names = ["LinQuant24","LogQuant16","RoundNearest16+Blosc","RoundNearest16+LZ4HC",
            "RoundNearest16+Zstd","RoundNearest16+Deflate","RoundNearest24+Blosc",
            "RoundNearest24+LZ4HC","RoundNearest24+Zstd","RoundNearest24+Deflate"]

# markers = ["s","s","d","d","d","d","o","o","o","o"]
markers = ["p","^","o","+","x","4","h","v"]

fig,ax = subplots(1,1,figsize=(8,6))
ax.set_xscale("log")
ax2 = ax.twinx()

ax2.plot(theoerror[1],1,"w",markersize=15,marker="*",lw=0,markeredgewidth=1.4,
            markeredgecolor="k",alpha=1,label="Float32",zorder=10)

for j in 1:length(allvars)
    ax.scatter(theoerror[1],LZs[j],120,color="w",marker=markers[j],
            edgecolor="C3",alpha=0.7)
end

for j in [1,8,2,7,3,4,5,6]
    if j in [4,5,6]
        ax.scatter([-1,-1,-1],[-1,-1,-1],120,color="k",edgecolor="k",marker=markers[j],
            alpha=0.55,lw=1.5,label=allvarsLN[j])
        ax.scatter(theoerror[1],LZs[j],120,color="C3",marker=markers[j],alpha=0.7)
    else
        ax.scatter([-1,-1,-1],[-1,-1,-1],120,color="w",edgecolor="k",marker=markers[j],
            alpha=0.55,lw=1.5,label=allvarsLN[j])
        ax.scatter(theoerror[1],LZs[j],120,color="w",marker=markers[j],
            edgecolor="C3",alpha=0.7)
    end
    ax.scatter(EVs[j,:,5],EVs[j,:,1],120,colours,marker=markers[j],
        alpha=0.55)
end

for i in [1,2,6,3,4,5,10,7,8,9]
    ax2.scatter(-1,-1,120,colours[i],marker="s",alpha=0.7,label=names[i])
end

ax2.plot(theoerror,theosize,"k",label="Theoretical rounding",zorder=-1)
ax.set_xlim(theoerror[1]/2,20)
ax.set_ylim(0.7,10)
ax2.set_ylim(0.7,10)
ax2.set_yticks([])
ax.set_yticks([1,2,4,6,8,10])

# ax.set_ylim(0,110)
ax.set_xlabel("max decimal error")
ax.set_ylabel("compression factor")
ax.legend(loc=2,ncol=1,fontsize=9,scatterpoints=1,title="Variable")
ax2.legend(loc=1,ncol=1,fontsize=9,scatterpoints=1,title="Compression method")
ax.set_title("CAMS data compression with various lossy/lossless methods",loc="left")

# ax.text(EV[1,3],EV[1,1],"median  ",rotation=90,ha="center",va="top",fontsize=9)
# ax.text(EV[1,4],EV[1,1],"90%  ",rotation=90,ha="center",va="top",fontsize=9)
# ax.text(EV[1,5],EV[1,1],"max  ",rotation=90,ha="center",va="top",fontsize=9)

tight_layout()
