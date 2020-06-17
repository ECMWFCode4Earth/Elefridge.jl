using NetCDF, LinearAlgebra
using Statistics, StatsBase
using PyCall
xr = pyimport("xarray")
#pygui(:qt)
using PyPlot
using MultivariateStats
# eofs = pyimport("eofs")

path = "/Users/milan/Downloads/cams/ieee"
filelist = filter(x->endswith(x,"no2.grib"),readdir(path))
gribfile = xr.open_dataset(joinpath(path,filelist[1]),engine="cfgrib")
no2 = gribfile.no2.values
level = 1:size(no2)[1]

no2m = mean(no2,dims=2)[:,1]
no2mi = median(no2,dims=2)[:,1]
no2min = minimum(no2,dims=2)[:,1]
no2max = maximum(no2,dims=2)[:,1]
no2p10 = [percentile(vec(no2[i,:]),10) for i in level]
no2p90 = [percentile(vec(no2[i,:]),90) for i in level]

## FIND FIRST SMALLER
function ispower2(x::Integer)
    while true
        if x == 2
            return true
        elseif x % 2 == 1
            return false
        else
            x = x ÷ 2
        end
    end
    return false
end


## FIND FIRST SMALLER
function findFirstSmaller(x::Float64,v::Array{Float64,1})

    l = length(v)-1
    @boundscheck ispower2(l) || throw(BoundsError())
    n = Int(log2(l))-1
    idx = l ÷ 2

    if x > v[end-1]     # round to max
        return l
    end

    # binary tree search
    for i in 1:n
        if x < v[idx]
            idx -= 2^(n-i)
        else
            idx += 2^(n-i)
        end
    end

    # split off the i = n+1 case
    if x >= v[idx]
        idx += 1
    end

    return idx-1
end

## LIN PACKING
Δlin24 = (Float64(maximum(no2max)) - Float64(minimum(no2min)))/(2^24-1)
Δlin16 = (Float64(maximum(no2max)) - Float64(minimum(no2min)))/(2^16-1)

lin24 = Array(minimum(no2min):Δlin24:maximum(no2max))
lin16 = Array(minimum(no2min):Δlin16:maximum(no2max))

lin24bounds = [(lin24[i+1]+lin24[i])/2 for (i,_) in enumerate(lin24[1:end-1])]
lin24bounds = vcat(minimum(no2min),lin24bounds,maximum(no2max))

lin16bounds = [(lin16[i+1]+lin16[i])/2 for (i,_) in enumerate(lin16[1:end-1])]
lin16bounds = vcat(minimum(no2min),lin16bounds,maximum(no2max))

no2lin24 = similar(no2)
no2lin16 = similar(no2)

for i in eachindex(no2)
    no2lin24[i] = Float32(lin24[findFirstSmaller(Float64(no2[i]),lin24bounds)])
    no2lin16[i] = Float32(lin16[findFirstSmaller(Float64(no2[i]),lin16bounds)])
end

## stats
no24m = mean(no2lin24,dims=2)[:,1]
no24mi = median(no2lin24,dims=2)[:,1]
no24min = minimum(no2lin24,dims=2)[:,1]
no24max = maximum(no2lin24,dims=2)[:,1]
no24p10 = [percentile(vec(no2lin24[i,:]),10) for i in level]
no24p90 = [percentile(vec(no2lin24[i,:]),90) for i in level]

no16m = mean(no2lin16,dims=2)[:,1]
no16mi = median(no2lin16,dims=2)[:,1]
no16min = minimum(no2lin16,dims=2)[:,1]
no16max = maximum(no2lin16,dims=2)[:,1]
no16p10 = [percentile(vec(no2lin16[i,:]),10) for i in level]
no16p90 = [percentile(vec(no2lin16[i,:]),90) for i in level]

## error 24bit
de24 = abs.(log10.(no2 ./ no2lin24))

de24m = mean(de24,dims=2)[:,1]
de24mi = median(de24,dims=2)[:,1]
de24min = minimum(de24,dims=2)[:,1]
de24max = maximum(de24,dims=2)[:,1]
de24p10 = [percentile(vec(de24[i,:]),10) for i in level]
de24p90 = [percentile(vec(de24[i,:]),90) for i in level]

de16 = abs.(log10.(no2 ./ no2lin16))

de16m = mean(de16,dims=2)[:,1]
de16mi = median(de16,dims=2)[:,1]
de16min = minimum(de16,dims=2)[:,1]
de16max = maximum(de16,dims=2)[:,1]
de16p10 = [percentile(vec(de16[i,:]),10) for i in level]
de16p90 = [percentile(vec(de16[i,:]),90) for i in level]

## plot
fig,(ax,ax2,ax3) = subplots(1,3,figsize=(10,6),sharey=true)
ax.invert_yaxis()

ax.semilogx(no24m,level,"C0",lw=3,label="mean")
ax.plot(no24mi,level,"C0",lw=1.5,label="median")
ax.fill_betweenx(level,no24min,no24max,color="C0",alpha=0.2,label="Min-max range")
ax.fill_betweenx(level,no24p10,no24p90,color="C0",alpha=0.5,label="10-90% range")

ax2.semilogx(no16m,level,"C0",lw=3,label="mean")
ax2.plot(no16mi,level,"C0",lw=1.5,label="median")
ax2.fill_betweenx(level,no16min,no16max,color="C0",alpha=0.2,label="Min-max range")
ax2.fill_betweenx(level,no16p10,no16p90,color="C0",alpha=0.5,label="10-90% range")

for i in level
    # 24-bit
    xvals = no2[i,(de24[i,:] .>= 1f-2) .* (de24[i,:] .< 1f-1)]
    ax.scatter(xvals,i*ones(length(xvals)),10,"yellow",alpha=0.05)

    xvals = no2[i,(de24[i,:] .>= 1f-1) .* (de24[i,:] .< 1f0)]
    ax.scatter(xvals,i*ones(length(xvals)),10,"C1",alpha=0.1)

    # 16-bit
    xvals = no2[i,(de16[i,:] .>= 2f-1) .* (de16[i,:] .< 0.5)]
    ax2.scatter(xvals,i*ones(length(xvals)),10,"yellow",alpha=0.05)

    xvals = no2[i,(de16[i,:] .>= 0.5) .* (de16[i,:] .< 1f0)]
    ax2.scatter(xvals,i*ones(length(xvals)),10,"C1",alpha=0.1)
end

ax.scatter([1,2,3],[0,0,0],10,"yellow",alpha=0.5,label="error > 1%")
ax.scatter([1,2,3],[0,0,0],10,"C1",alpha=0.5,label="error > 10%")

ax2.scatter([1,2,3],[0,0,0],10,"yellow",alpha=0.5,label="error > 20%")
ax2.scatter([1,2,3],[0,0,0],10,"C1",alpha=0.5,label="error > 50%")

ax3.semilogx(de24m,level,"C2",lw=3,label="24bit")
ax3.plot(de24mi,level,"C2",lw=1.5)
ax3.fill_betweenx(level,de24min,de24max,color="C2",alpha=0.2)
ax3.fill_betweenx(level,de24p10,de24p90,color="C2",alpha=0.5)

ax3.semilogx(de16m,level,"C1",lw=3,label="16bit")
ax3.plot(de16mi,level,"C1",lw=1.5)
ax3.fill_betweenx(level,de16min,de16max,color="C1",alpha=0.2)
ax3.fill_betweenx(level,de16p10,de16p90,color="C1",alpha=0.5)

ax.set_title(L"NO$_2$: 24-bit linear packing")
ax2.set_title(L"NO$_2$: 16-bit linear packing")
ax3.set_title("L1 log-error")
ax.set_ylabel("model level")
ax.set_xlabel("mixing ratio kg/kg")
ax2.set_xlabel("mixing ratio kg/kg")
ax.set_ylim(137,1)
ax.legend(loc=2,scatterpoints=3)
ax2.legend(loc=2,scatterpoints=3)
ax3.legend(loc=1)

ax.set_xlim(1e-14,1e-6)
ax2.set_xlim(1e-14,1e-6)

#ax2.set_xticks(lin16[1:20])

tight_layout()
