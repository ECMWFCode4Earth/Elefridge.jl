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

## log PACKING
Δlog24 = (log10(Float64(maximum(no2max))) - log10.(Float64(minimum(no2min))))/(2^24-1)
Δlog16 = (log10(Float64(maximum(no2max))) - log10.(Float64(minimum(no2min))))/(2^16-1)

log24 = 10.0 .^ Array(log10(minimum(no2min)):Δlog24:log10(maximum(no2max)))
log16 = 10.0 .^ Array(log10(minimum(no2min)):Δlog16:log10(maximum(no2max)))

log24bounds = [(log24[i+1]+log24[i])/2 for (i,_) in enumerate(log24[1:end-1])]
log24bounds = vcat(minimum(no2min),log24bounds,maximum(no2max))

log16bounds = [(log16[i+1]+log16[i])/2 for (i,_) in enumerate(log16[1:end-1])]
log16bounds = vcat(minimum(no2min),log16bounds,maximum(no2max))

no2log24 = similar(no2)
no2log16 = similar(no2)

for i in eachindex(no2)
    no2log24[i] = Float32(log24[findFirstSmaller(Float64(no2[i]),log24bounds)])
    no2log16[i] = Float32(log16[findFirstSmaller(Float64(no2[i]),log16bounds)])
end

## stats
no24m = mean(no2log24,dims=2)[:,1]
no24mi = median(no2log24,dims=2)[:,1]
no24min = minimum(no2log24,dims=2)[:,1]
no24max = maximum(no2log24,dims=2)[:,1]
no24p10 = [percentile(vec(no2log24[i,:]),10) for i in level]
no24p90 = [percentile(vec(no2log24[i,:]),90) for i in level]

no16m = mean(no2log16,dims=2)[:,1]
no16mi = median(no2log16,dims=2)[:,1]
no16min = minimum(no2log16,dims=2)[:,1]
no16max = maximum(no2log16,dims=2)[:,1]
no16p10 = [percentile(vec(no2log16[i,:]),10) for i in level]
no16p90 = [percentile(vec(no2log16[i,:]),90) for i in level]

## error 24bit
de24 = abs.(log10.(no2 ./ no2log24))

de24m = mean(de24,dims=2)[:,1]
de24mi = median(de24,dims=2)[:,1]
de24min = minimum(de24,dims=2)[:,1]
de24max = maximum(de24,dims=2)[:,1]
de24p10 = [percentile(vec(de24[i,:]),10) for i in level]
de24p90 = [percentile(vec(de24[i,:]),90) for i in level]

# and with 16
de16 = abs.(log10.(no2 ./ no2log16))

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

# for i in level
#     # 24-bit
#     xvals = no2[i,(de24[i,:] .>= 3f-7) .* (de24[i,:] .< 1f0)]
#     ax.scatter(xvals,i*ones(length(xvals)),10,"yellow",alpha=0.05)
#
#     # 16-bit
#     xvals = no2[i,(de16[i,:] .>= 5.45f-5) .* (de16[i,:] .< 1f0)]
#     ax2.scatter(xvals,i*ones(length(xvals)),10,"yellow",alpha=0.05)
# end

# ax.scatter([1,2,3],[0,0,0],10,"yellow",alpha=0.5,label="error > 5e-7")
# ax2.scatter([1,2,3],[0,0,0],10,"yellow",alpha=0.5,label="error > 0.006%")

ax3.semilogx(de24m,level,"C2",lw=3,label="24bit")
ax3.plot(de24mi,level,"C2",lw=1.5)
ax3.fill_betweenx(level,de24min,de24max,color="C2",alpha=0.2)
ax3.fill_betweenx(level,de24p10,de24p90,color="C2",alpha=0.5)

ax3.semilogx(de16m,level,"C1",lw=3,label="16bit")
ax3.plot(de16mi,level,"C1",lw=1.5)
ax3.fill_betweenx(level,de16min,de16max,color="C1",alpha=0.2)
ax3.fill_betweenx(level,de16p10,de16p90,color="C1",alpha=0.5)

ax.set_title(L"NO$_2$: 24-bit log packing")
ax2.set_title(L"NO$_2$: 16-bit log packing")
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

ax3.set_xticks([1e-8,1e-7,1e-6,1e-5,1e-4])

tight_layout()
