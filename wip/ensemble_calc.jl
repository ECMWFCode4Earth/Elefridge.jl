using JLD2, FileIO, NetCDF
using Elefridge

path = "/network/aopp/chaos/pred/kloewer/esowc/"
member1files = filter(x->endswith(x,".nc"),readdir(joinpath(path,"member1")))
steps = [parse(Int,split(split(file,"step")[2],".")[1]) for file in member1files]
sort!(steps)

# parameter
Nens = 25
Nbits = 32
Nlon = 1800
Nlat = 901
Nvert = 91
Ntsteps = length(steps)

# preallocate/load bitwise information content
BI = zeros(Float64,Ntsteps,Nbits)
try
    BI[:,:] = load(joinpath(path,"analysis/bitinformation_$(Nens)members.jld2"),"BI")
    println("Result array loaded from file.")
catch
    println("Result array preallocated.")
end

temp = fill(0f0,Nens,Nlon,Nlat,Nvert)

for t in 1:Ntsteps
    println("Time step $(steps[t])h")
    
    # read the data
    print("Ensemble member ")
    for ie in 1:Nens
        print("$ie,")
        ncfile = NetCDF.open(joinpath(path,"member$ie","ensemble.t.member$ie.step$(steps[t]).ll.nc"))
        temp[ie,:,:,:] = ncfile.vars["t"][1:Nlon,1:Nlat,:]
    end
    println("All loaded.")
    
    # calculate bitwise information
    BI[t,:] = bitinformation(temp)
    println("Information calculated.")
    save(joinpath(path,"analysis/bitinformation_$(Nens)members.jld2"),"BI",BI)
    println("---")
end