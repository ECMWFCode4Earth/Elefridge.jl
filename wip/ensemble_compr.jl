using JLD2, FileIO, NetCDF
using Elefridge, ZfpCompression
using TranscodingStreams,CodecZstd

ZstdCompressorL3 = ZstdCompressor(level=3)
TranscodingStreams.initialize(ZstdCompressorL3)

ZstdCompressorL10 = ZstdCompressor(level=10)
TranscodingStreams.initialize(ZstdCompressorL10)

ZstdCompressorL22 = ZstdCompressor(level=22)
TranscodingStreams.initialize(ZstdCompressorL22)

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
chunk_size = 8
nchunks = Nvert ÷ chunk_size
keepbits = [8,7,6,5]

# preallocate
Xzfp,Xzstd = try load(joinpath(path,"analysis/ensemble_compression_$(Nens)members_s5-8.jld2"),"Xzfp","Xzstd") catch; zeros(Float64,Ntsteps,length(keepbits)),zeros(Float64,Ntsteps,length(keepbits)) end
    
for t in 13:Ntsteps
    println("Time step $(steps[t])h")
    temp = fill(0f0,Nens,Nlon,Nlat,Nvert)
   
    print("Loading ensemble member ")
    # read the data
    for ie in 1:Nens
        print("$ie,")
        ncfile = NetCDF.open(joinpath(path,"member$ie","ensemble.t.member$ie.step$(steps[t]).ll.nc"))
        temp[ie,:,:,:] = ncfile.vars["t"][1:Nlon,1:Nlat,1:Nvert]
    end
    println("\nCompressing...")
    
    # ZFP COMPRESSION
    for (ik,keepbit) in enumerate(keepbits)
        
        precision = keepbit+7   # For 4D compression
    
        # ZFP: compress in chunks
        compressed_size = 0
    
        for ichunk in 1:nchunks-1    # compress the last one after loop
            print("$ichunk,")
            start = (ichunk-1)*chunk_size+1
            stop = start+chunk_size-1
            temp_chunk = temp[:,:,:,start:stop]
            compressed_size += sizeof(zfp_compress(temp_chunk,precision=precision))
        end
    
        # compress last chunk
        start = (nchunks-1)*chunk_size+1
        temp_chunk = temp[:,:,:,start:end]
        compressed_size += sizeof(zfp_compress(temp_chunk,precision=precision))
    
        Xzfp[t,ik] = 2*sizeof(temp)/compressed_size
        println("$(keepbit)bits, CF ZFP: $(Xzfp[t,ik])")
        save(joinpath(path,"analysis/ensemble_compression_$(Nens)members_s5-8.jld2"),"Xzfp",Xzfp,"Xzstd",Xzstd)
    end
    
    # ZSTD COMPRESSION
    for (ik,keepbit) in enumerate(keepbits)
        round!(temp,keepbit)
        Xr8 = copy(reinterpret(UInt8,vec(temp)))
        Xzstd[t,ik] = 2*sizeof(Xr8)/sizeof(transcode(ZstdCompressorL10,Xr8))
        println("$(keepbit)bits, CF ZSTD: $(Xzstd[t,ik])")
        save(joinpath(path,"analysis/ensemble_compression_$(Nens)members_s5-8.jld2"),"Xzfp",Xzfp,"Xzstd",Xzstd)
    end
    
    println("Compression done.")
    println("---")
end
