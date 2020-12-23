using PyCall
using JLD2
using FileIO
using NetCDF
using Elefridge

steps = "1-90"
# steps = "90-144"
# steps = "150-360"
path1 = "/data1/ens"
path2 = "/data2/ens"
all_files1 = [joinpath(path1,f) for f in filter(x->endswith(x,"steps$steps.ll.nc"),readdir(path1))]
all_files2 = [joinpath(path2,f) for f in filter(x->endswith(x,"steps$steps.ll.nc"),readdir(path2))]

all_files = vcat(all_files1,all_files2)

# sort members
all_files_n = sortperm([parse(Int,split(split(f,"member")[2],".steps")[1])
                        for f in all_files])
all_files = all_files[all_files_n]
all_files

# parameter
Nens = length(all_files)
Ntsteps = [90,18,36][steps .== ["1-90","90-144","150-360"]][1]
Nbits = 32
Nlon = 1800
Nlat = 451
Nvert = 11

# preallocate bitwise information content
BI = try load("/home/milan/analysis/bitinformation_$(Nens)members_$steps.jld2","BI") catch; zeros(Float64,Ntsteps,Nbits) end

temp = fill(0f0,Nens,Nlon,Nlat,Nvert)

for t in 1:Ntsteps
    println("Time step $t")
    temp[:] .= 0f0
    
    # read the data
    print("Ensemble member ")
    for ie in 1:Nens
        print("$ie,")
        ncfile = NetCDF.open(all_files[ie])
        temp[ie,:,:,:] = ncfile.vars["t"][:,451:end,:,t]
    end
    println("All loaded.")
    
    # calculate bitwise information
    BI[t,:] = bitinformation(temp)
    println("Information calculated.")
    save("/home/milan/analysis/bitinformation_$(Nens)members_$steps.jld2","BI",BI)
    println("---")
end