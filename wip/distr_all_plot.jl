using PyPlot
using JLD
using Printf

@load "/Users/milan/cams/entropy/distribution_all_gridded.jld"

n = length(S)

## sort and group
aero = Array(1:15)
ozone = [35,44,45]
methane = [25,26,41]
dynamics = [33,34,54,55]
clouds = [27,28,31,32,51]
hydro = vcat(Array(36:40),46)
nitro = [42,43]
oopp = Array(47:50)
ces = vcat(Array(17:21),[23,24])

groups = [aero,ozone,methane,dynamics,clouds,hydro,nitro,oopp,ces]
grouped = vcat(aero,ozone,methane,dynamics,clouds,hydro,nitro,oopp,ces)

n_in_groups = n-length(grouped)
others = fill(0,n_in_groups)
j = 1
for i in 1:55
    if i in grouped
        0
    else
        global others[j] = i
        global j+= 1
    end
end

ax1names = fill("",55)
ax2names = fill("",55)

groups1 = [aero,nitro,hydro,ces]
groups2 = [dynamics,methane,clouds,oopp,ozone,others]

j = 36
for group in groups1
    for i in group
        global ax1names[j] = varnames[i]
        global j -= 1
    end
    global j -= 2
end

j = 35
for group in groups2
    for i in group
        global ax2names[j] = varnames[i]
        global j -= 1
    end
    global j -= 2
end


## plotting
pygui(true)
fig,(ax1,ax2) = subplots(1,2,figsize=(10,8))
ax1.set_xscale("log")
ax2.set_xscale("log")

ax1y2 = ax1.twinx()
ax2y2 = ax2.twinx()

groups1colors = ["C4","indianred","C9","C8"]
groups2colors = ["C2","C1","royalblue","C5","C6","darkslategray"]
ioffset = 35

for (ig,group) in enumerate(groups1)

    color = groups1colors[ig]

    for i in group
        global yoffset = 0.01*ioffset

        x = (S[i].edges[1][4:end-3] + S[i].edges[1][5:end-2])/2
        y = yoffset .+ S[i].weights[4:end-1]

        # extend line to zero
        x = vcat(0,S[i].edges[1][4],x,S[i].edges[1][end-2:end-1])
        y = vcat(yoffset,yoffset,y,S[i].weights[end-1]+yoffset)

        l = length(y)

        ax1.plot(x,y,color)
        ax1.fill_between(x,y,fill(yoffset,l),color=color,alpha=0.3)
        global ioffset -= 1
    end
    global ioffset -= 2
end

ioffset = 34

for (ig,group) in enumerate(groups2)

    color = groups2colors[ig]

    for i in group
        global yoffset = 0.01*ioffset

        x = (S[i].edges[1][4:end-3] + S[i].edges[1][5:end-2])/2
        y = yoffset .+ S[i].weights[4:end-1]

        # extend line to zero
        x = vcat(0,S[i].edges[1][4],x,S[i].edges[1][end-2:end-1])
        y = vcat(yoffset,yoffset,y,S[i].weights[end-1]+yoffset)

        l = length(y)

        ax2.plot(x,y,color)
        ax2.fill_between(x,y,fill(yoffset,l),color=color,alpha=0.3)
        global ioffset -= 1
    end
    global ioffset -= 2
end

ax1.set_xlim(1e-26,1e-6)
ax1.set_ylim(-0.01,0.4)
ax1y2.set_ylim(-0.01,0.4)

ax2.set_xlim(1e-26,1e3)
ax2.set_ylim(-0.01,0.4)
ax2y2.set_ylim(-0.01,0.4)

ytik = Array(0:0.05:0.4)
ytikm = Array(0:0.01:0.4)

ax1.set_yticks(ytik)
ax1.set_yticklabels([@sprintf("%i%%",y*100) for y in ytik])
ax1.set_yticks(ytikm,minor=true)

ax2.set_yticks(ytik)
ax2.set_yticklabels([])
ax2.set_yticks(ytikm,minor=true)

ax1y2.set_yticks(ytikm)
ax1y2.set_yticklabels(ax1names)
ax2y2.set_yticks(ytikm)
ax2y2.set_yticklabels(ax2names)

ax1.set_xlabel("value")
ax2.set_xlabel("value")
ax1.set_ylabel("frequency")

tight_layout()
