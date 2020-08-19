function bittranspose(A::Array{T,1}) where {T<:Unsigned}

    nbits = sizeof(eltype(A))*8
    N = length(A)

    B = fill(zero(T),N)                 	# preallocate
    n_in_block = 2^10                  		# values per block
    nblocks = (N-1) ÷ n_in_block + 1		# number of blocks

    for nb in 1:nblocks

        lastindex = min(nb*n_in_block,N)

        Ablock = view(A,(nb-1)*n_in_block+1:lastindex)
        Bblock = view(B,(nb-1)*n_in_block+1:lastindex)

        i = 0   # counter for transposed bits
        for bi in 1:nbits
            # mask to extract bit bi in each element of A
            mask = one(T) << (nbits-bi)

			# walk through elements in A first, then change bit-position
			# this comes with disadvantage that A has to be read nbits-times
			# but allows for more natural indexing, as all the sign bits are
			# read first, and so on.

            for (ia,a) in enumerate(Ablock)
                # mask non-bi bits and
                # (1) shift by nbits-bi >> to the right, either 0x0 or 0x1
                # (2) shift by (nbits-1) - (i % nbits) << to the left
                # combined this is: >> ((i % nbits)-bi+1)
                bit = (a & mask) >> ((i % nbits)-bi+1)

                # the first nbits elements go into same b in B, then next b
                @inbounds Bblock[(i ÷ nbits) + 1] |= bit
                i += 1
            end
        end
    end

    return B
end

# function bittranspose(::Type{T},A::Array{T2,1}) where {T<:Unsigned,T2<:AbstractFloat}
#     B = reinterpret.(T,A)
#     At = bittranspose(A_asUInt)
#     At = unsafe_wrap(Array, Ptr{T2}(pointer(At)), sizeof(At))
#     return At
# end

bittranspose(A::Array{Float16,1}) = reinterpret.(Float16,bittranspose(reinterpret.(UInt16,A)))
bittranspose(A::Array{Float32,1}) = reinterpret.(Float32,bittranspose(reinterpret.(UInt32,A)))
bittranspose(A::Array{Float64,1}) = reinterpret.(Float64,bittranspose(reinterpret.(UInt64,A)))

bittranspose(A::AbstractArray) where {T,N} = reshape(bittranspose(vec(A)),size(A))

function bitbacktranspose(A::Array{T,1}) where {T<:Unsigned}

    nbits = sizeof(eltype(A))*8
    N = length(A)

    B = fill(zero(T),N)                 # preallocate
    n_in_block = 2^10                  	# values per block
    nblocks = (N-1) ÷ n_in_block + 1	# number of blocks

    for nb in 1:nblocks

        lastindex = min(nb*n_in_block,N)

        Ablock = view(A,(nb-1)*n_in_block+1:lastindex)
        Bblock = view(B,(nb-1)*n_in_block+1:lastindex)

		nelements = length(Bblock)	# = n_in_block except for the last
									# block where usually smaller

        i = 0   # counter for transposed bits
		for (ia,a) in enumerate(Ablock)
        	for bi in 1:nbits
            	# mask to extract bit bi in each element of A
            	mask = one(T) << (nbits-bi)

                # mask non-bi bits and
                # (1) shift by nbits-bi >> to the right, either 0x0 or 0x1
                # (2) shift by (nbits-1) - (i ÷ nblockfloat) << to the left
                # combined this is: >> ((i % nbits)-bi+1)
                bit = (a & mask) >> ((i ÷ nelements)-bi+1)

                # the first nbits elements go into same b in B, then next b
                @inbounds Bblock[(i % nelements) + 1] |= bit
                i += 1
            end
        end
    end

    return B
end

bitbacktranspose(A::Array{Float16,1}) = reinterpret.(Float16,bitbacktranspose(reinterpret.(UInt16,A)))
bitbacktranspose(A::Array{Float32,1}) = reinterpret.(Float32,bitbacktranspose(reinterpret.(UInt32,A)))
bitbacktranspose(A::Array{Float64,1}) = reinterpret.(Float64,bitbacktranspose(reinterpret.(UInt64,A)))

bitbacktranspose(A::AbstractArray) where {T,N} = reshape(bitbacktranspose(vec(A)),size(A))

# function Base.BitMatrix(A::Array{T,1}) where T
# 	isbitstype(eltype(A)) || error("Only bitstype for elements of A allowed.
# 									$(eltype(A)) provided")
# 	height = 8 * sizeof(eltype(A))
# 	dims = (height, length(A))
# 	data_elems = cld(sizeof(A), 8)
# 	bitarr = Base.BitMatrix(undef, dims)
# 	bitarr.chunks = unsafe_wrap(Array, Ptr{UInt}(pointer(A)), (data_elems,))
# 	return bitarr
# end
#
# function foo(x::Array)
# 	isbitstype(eltype(x)) || error("Bad!")
# 	height = 8 * sizeof(eltype(x))
# 	dims = (height, length(x))
# 	data_elems = cld(sizeof(x), 8)
# 	bitarr = Base.BitMatrix(undef, dims)
# 	bitarr.chunks = unsafe_wrap(Array, Ptr{UInt}(pointer(x)), (data_elems,))
# 	return bitarr
# end
