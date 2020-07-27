function bitentropy(A::AbstractArray,base::Real=2)
    nbits = sizeof(eltype(A))*8
    T = whichUInt(nbits)
    nbits == 8 && return bitentropy(T,A,base)
    nbits == 16 && return bitentropy(T,A,base)
    nbits == 24 && return bitentropy(T,A,base)
    nbits == 32 && return bitentropy(T,A,base)
    nbits == 40 && return bitentropy(T,A,base)
    nbits == 48 && return bitentropy(T,A,base)
    nbits == 56 && return bitentropy(T,A,base)
    nbits == 64 && return bitentropy(T,A,base)
    throw(error("Only element types with 8,16,24,32,40,48,56,64 bits supported.
                $nbits-bit $(eltype(A)) provided."))
end

function bitentropy(::Type{T},A::AbstractArray,base::Real=2) where T

    # reinterpret to UInt then sort to avoid allocating a histogram
    Av = sort(reinterpret.(T,vec(A)))
    n = length(Av)

    E = 0.0     # entropy
    m = 1       # counter

    for i in 1:n-1
        @inbounds if Av[i] == Av[i+1]     # check whether next entry belongs to the same bin in histogram
            m += 1
        else
            p = m/n
            E -= p*log(p)
            m = 1
        end
    end

    p = m/n         # complete for last bin
    E -= p*log(p)

    # convert to given base, 2 i.e. [bit] by default
    E /= log(base)

    return E
end
