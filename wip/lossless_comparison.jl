using PyPlot
using JLD

path = "/Users/milan/cams"
allvars = ["no2","go3","so2","aermr04","aermr05","aermr06","ch4","co","co2","no"]
vari = "so2"

# volume and error norms (columns):
# vol, mean, median, 90% and max  of decimal error
# mean, median, 90% and max of abs error
# total 9

# compression methods (rows):
# LinQuant24, LogQuant16, RoundNearest16/24 + Blosc/Deflate/Zstd/LZ4
# total 10

D = load("cams/error/$vari.jld")
EV = D["EV"]
lz4hcc = D["lz4"]

# theoretical
theosize = [32/i for i in 32:-1:5]
x0 = 1f0
theoerror = [log10(nextfloat(x0,2^i)) for i in 0:27]

## plotting
pygui(true)
n = size(EV)[1]

cm1 = ColorMap("viridis")
cm2 = ColorMap("inferno")

function RGB(cv::Tuple)
    r = UInt8(round(cv[1]*255))
    g = UInt8(round(cv[2]*255))
    b = UInt8(round(cv[3]*255))
    return "#"*repr(r)[3:4]*repr(g)[3:4]*repr(b)[3:4]
end

colours = ["0.5","0.3",RGB(cm2(0.9)),RGB(cm2(0.8)),RGB(cm2(0.6)),RGB(cm2(0.4)),
                    RGB(cm1(0.8)),RGB(cm1(0.6)),RGB(cm1(0.4)),RGB(cm1(0.2))]

names = ["LinQuant24","LogQuant16","RoundNearest16+Blosc","RoundNearest16+LZ4HC",
            "RoundNearest16+Zstd","RoundNearest16+Deflate","RoundNearest24+Blosc",
            "RoundNearest24+LZ4HC","RoundNearest24+Zstd","RoundNearest24+Deflate"]

markers = ["s","s","d","d","d","d","o","o","o","o"]

fig,ax = subplots(1,1,figsize=(8,6))
ax.set_xscale("log")

ax.plot(theoerror[1],1,"w",markersize=15,marker="*",lw=0,markeredgewidth=1.4,
            markeredgecolor="k",alpha=1,label="Float32",zorder=10)
ax.plot(theoerror[1],lz4hcc,"w",markersize=10,marker="v",lw=0,markeredgewidth=1.4,
            markeredgecolor="k",alpha=0.7,label="Float32+LZ4HC",zorder=10)

for i in [1,2,3,4,5,6,7,8,9,10]
    # ax.plot(EV[i,[3,5]],[EV[i,1],EV[i,1]],"k",lw=5.5,zorder=-1)
    ax.plot(EV[i,3:4],[EV[i,1],EV[i,1]],colours[i],
            marker=markers[i],markevery=2,
            markeredgecolor="k",markersize=10,
            lw=5,label=names[i])
    ax.plot(EV[i,4:5],[EV[i,1],EV[i,1]],colours[i],lw=2)
end

ax.plot(theoerror,theosize,"k",label="Theoretical rounding",zorder=-1)
ax.set_xlim(theoerror[1]/3,100)
ax.set_ylim(0.3,1+ceil(maximum(EV[:,1])))
# ax.set_ylim(0,110)
ax.set_xlabel("decimal error")
ax.set_ylabel("compression factor")
ax.legend(loc=1,ncol=1,fontsize=9)
ax.set_title("NO compression",loc="left")

ax.text(EV[1,3],EV[1,1],"median  ",rotation=90,ha="center",va="top",fontsize=9)
ax.text(EV[1,4],EV[1,1],"90%  ",rotation=90,ha="center",va="top",fontsize=9)
ax.text(EV[1,5],EV[1,1],"max  ",rotation=90,ha="center",va="top",fontsize=9)

tight_layout()
