using BitInformation
using JLD2, NetCDF

path = "/network/aopp/chaos/pred/kloewer/esowc/cams/"
filelist = filter(x->endswith(x,"_v3.nc"),readdir(path))

n = length(filelist)
varnames = fill("",n)
nbits = 32

IC = fill(0.0,n,nbits)

for (i,file) in enumerate(filelist)
    varname = split(split(file,"cams_")[end],"_2019")[1]
    varnames[i] = varname
    print("$varname")
    ncfile = NetCDF.open(joinpath(path,file))
    
    # find the variable name by size
    var = [var for var in ncfile.vars if prod(size(var[2])) == 900*451*137][1][1]
    X = ncfile.vars[var][:,6:end-5,:]  # exclude north/southpole
    
    # convert biased exponent to signed exponent
    BitInformation.signed_exponent!(X)

    IC[i,:] = bitinformation(X,:all_dimensions)
    print("X, ")    

    @save joinpath(path,"..","analysis/bitinformation_all_filt.jld2") varnames IC
end

