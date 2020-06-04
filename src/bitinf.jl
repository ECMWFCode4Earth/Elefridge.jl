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

# Float32 information entropy
bitmin = Int(reinterpret(UInt32,minimum(no2)))
bitmax = Int(reinterpret(UInt32,maximum(no2)))
diff = bitmax-bitmin
B = Array{UInt16,1}(undef,diff+1)

for x in no2
    B[reinterpret(UInt32,x)-bitmin+1] += 0x1
end

# throw away all zero entries (don't contribute to entropy)
B = B[B .!= 0x0000]
p = B / sum(B)

H_F32 = entropy(p,2)

# Lin24 information entropy
L = LinQuant24Array(no2)
B = zeros(UInt16,2^24)

for x in L
    B[x+1] += 0x1
end

B = B[B .!= 0x0000]
p = B / sum(B)

H_Lin24 = entropy(p,2)

# Log16 information entropy
L = LogQuant16Array(no2)
B_log = zeros(UInt16,2^16)

for x in L
    B_log[x+1] += 0x1
end

p = B_log / sum(B_log)
H_Log16 = entropy(p,2)
