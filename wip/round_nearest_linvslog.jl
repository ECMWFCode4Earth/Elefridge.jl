using Statistics
using StatsBase
using Elefridge
using PyCall
xr = pyimport("xarray")
using JLD

path = "/Users/milan/cams/gridded"
filelist = filter(x->endswith(x,".grib"),readdir(path))

nvars = length(filelist)
nmethods = 4        # Lin24,Lin16,Log16lin,Log16log
nstats = 7          # min,10%,25%,50%,75%,90%,max

meanerror = fill(0f0,nmethods,nvars)
linerror = fill(0f0,nmethods,nvars,nstats)
logerror = fill(0f0,nmethods,nvars,nstats)
allvars = Array{String,1}(undef,nvars)

for i in 1:nvars

    file = filelist[i]
    allvars[i] = split(split(file,"_")[end],".")[1]

    println("---")
    println("Compressing $(allvars[i])")

    D = xr.open_dataarray(joinpath(path,file),engine="cfgrib").data

    if any(D .<= 0f0)
        println("Negative/zero entries found.")
        D = D[D.>0f0]      # remove negative entries
    end

    if any(isnan.(D))
        println("NaNs found.")
        D = D[.~isnan(D)]   # remove NaNs
    end

    Dlin16 = Array(LinQuant16Array(D))
    Dlin24 = Array(LinQuant24Array(D))
    Dlog16log = Array(LogQuant16Array(D,:logspace))
    Dlog16lin = Array(LogQuant24Array(D,:logspace))

    # mean error
    meanerror[1,i] = abs(mean(D - Dlin16))
    meanerror[2,i] = abs(mean(D - Dlin24))
    meanerror[3,i] = abs(mean(D - Dlog16log))
    meanerror[4,i] = abs(mean(D - Dlog16lin))

    # absolute error
    for (j,Dr) in enumerate([Dlin16,Dlin24,Dlog16log,Dlog16lin])
        E = sort(vec(abs.(D-Dr)))
        linerror[j,i,1],linerror[j,i,7] = E[1],E[end]
        linerror[j,i,2] = quantile(E,0.1,sorted=true)
        linerror[j,i,3] = quantile(E,0.25,sorted=true)
        linerror[j,i,4] = quantile(E,0.5,sorted=true)
        linerror[j,i,5] = quantile(E,0.75,sorted=true)
        linerror[j,i,6] = quantile(E,0.9,sorted=true)
    end

    # decimal error
    for (j,Dr) in enumerate([Dlin16,Dlin24,Dlog16log,Dlog16lin])
        E = sort(vec(abs.(log10.(abs.(D./Dr)))))    # avoid non-finite entries with abs
        logerror[j,i,1],logerror[j,i,7] = E[1],E[end]
        logerror[j,i,2] = quantile(E,0.1,sorted=true)
        logerror[j,i,3] = quantile(E,0.25,sorted=true)
        logerror[j,i,4] = quantile(E,0.5,sorted=true)
        logerror[j,i,5] = quantile(E,0.75,sorted=true)
        logerror[j,i,6] = quantile(E,0.9,sorted=true)
    end
end

save("/Users/milan/cams/error/linvslog_all.jld",
        "allvars",allvars,
        "meanerror",meanerror,
        "linerror",linerror,
        "logerror",logerror)
