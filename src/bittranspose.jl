function bittranspose(A::Array{T,1}) where {T<:Unsigned}

    nbits = sizeof(eltype(A))*8
    N = length(A)

    B = fill(zero(T),N)                 # preallocate
    nblockfloat = 2^10                  # float32 values per block
    nblockbyte = nblockfloat*4          # blocksize in byte
    nblockbit = nblockbyte*8            # blocksize in bit

    nblocks = (N-1) ÷ nblockbyte + 1

    for nb in 1:nblocks

        lastindex = min(nb*nblockfloat,N)

        Ablock = view(A,(nb-1)*nblockfloat+1:lastindex)
        Bblock = view(B,(nb-1)*nblockfloat+1:lastindex)

        i = 0   # counter for transposed bits
        for bi in 1:nbits
            # mask to extract bit bi in each element of A
            mask = one(T) << (nbits-bi)

            for (ia,a) in enumerate(Ablock)
                ui = reinterpret(T,a)

                # mask non-bi bits and
                # (1) shift by nbits-bi >> to the right, either 0x0 or 0x1
                # (2) shift by (nbits-1) - (i % nbits) << to the left
                # combined this is: >> ((i % nbits)-bi+1)
                bit = (ui & mask) >> ((i % nbits)-bi+1)

                # the first nbits elements go into same b in B, then next b
                Bblock[(i ÷ nbits) + 1] |= bit
                i += 1
            end
        end
    end

    return B
end

bittranspose(A::Array{Float16,1}) = reinterpret.(Float16,bittranspose(reinterpret.(UInt16,A)))
bittranspose(A::Array{Float32,1}) = reinterpret.(Float32,bittranspose(reinterpret.(UInt32,A)))
bittranspose(A::Array{Float64,1}) = reinterpret.(Float64,bittranspose(reinterpret.(UInt64,A)))

bittranspose(A::AbstractArray) where {T,N} = reshape(bittranspose(vec(A)),size(A))
