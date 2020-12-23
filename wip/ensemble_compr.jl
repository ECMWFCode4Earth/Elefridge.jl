using PyCall
using JLD2
using FileIO
using NetCDF
using Elefridge
using ZfpCompression
using TranscodingStreams,CodecZstd

ZstdCompressorL3 = ZstdCompressor(level=3)
TranscodingStreams.initialize(ZstdCompressorL3)

ZstdCompressorL10 = ZstdCompressor(level=10)
TranscodingStreams.initialize(ZstdCompressorL10)

ZstdCompressorL22 = ZstdCompressor(level=22)
TranscodingStreams.initialize(ZstdCompressorL22)

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
Nlon = 1800
Nlat = 451
Nvert = 11
every = 5

# preallocate bitwise information content
Xzfp,Xzstd = try load("/home/milan/analysis/ensemble_compression_$(steps).jld2","Xzfp","Xzstd") catch; zeros(Float64,Ntsteps),zeros(Float64,Ntsteps) end
    
for t in 77:Ntsteps
    println("Time step $t")
    temp = fill(0f0,Nens,Nlon,Nlat,Nvert)
   
    print("Ensemble member ")
    # read the data
    for ie in 1:Nens
        print("$ie,")
        ncfile = NetCDF.open(all_files[ie])
        temp[ie,:,:,:] = ncfile.vars["t"][:,451:end,:,t]
    end
    println("All loaded.")
    
#     precision = t<5 ? 14 : 13
#     keepbits = t<5 ? 7 : 6
    precision = 14
    keepbits = 7
    
    # compress
    Xzfp[t] = 2*sizeof(temp)/sizeof(zfp_compress(temp,precision=precision))
    println("CF ZFP: $(Xzfp[t])")

    round!(temp,keepbits)
    Xr8 = copy(reinterpret(UInt8,vec(temp)))
    Xzstd[t] = 2*sizeof(Xr8)/sizeof(transcode(ZstdCompressorL10,Xr8))
    println("CF ZSTD: $(Xzstd[t])")
    
    println("Compression done.")
    save("/home/milan/analysis/ensemble_compression_$(steps).jld2","Xzfp",Xzfp,"Xzstd",Xzstd)
    println("---")
end