function bitentropy(n::Integer,A::AbstractArray)

    n == 8*sizeof(eltype(A)) || throw(error("Element type is $(eltype(A)) with
                            $(sizeof(eltype(A))) bits, incompatible with n=$n"))

    T = whichUInt(n)
    hist = zeros(UInt32,2^n)

    for a in A
        hist[reinterpret(T,a)+1] += one(UInt32)
    end

    # normalize
    hist /= sum(hist)
    return entropy(hist,2)
end

function bitentropy(A::AbstractArray)
    nbits = sizeof(eltype(A))*8
    nbits == 8 && return bitentropy(8,A)
    nbits == 16 && return bitentropy(16,A)
    nbits == 24 && return bitentropy(24,A)
    throw(error("Only element types with 8,16,24 bits supported.
                $nbits-bit $(eltype(A)) provided."))
end
