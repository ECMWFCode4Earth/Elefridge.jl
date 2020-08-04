function bitcount(A::Array{T},b::Int) where {T<:Unsigned}
    N = sizeof(T)*8             # number of bits in T
    @boundscheck b <= N || throw(BoundsError("Count bit $b for $T is invalid."))
    n = 0                       # counter
    shift = N-b                 # shift desired bit b
    mask = one(T) << shift
    for a in A                      # mask everything but b and shift
        n += (a & mask) >> shift    # to have either 0x00 or 0x01
    end
    return n
end

# function bitpaircount(A::Array{T},b::Int) where {T<:Unsigned}
#     N = sizeof(T)*8             # number of bits in T
#     @boundscheck b <= N || throw(BoundsError("Count bit $b for $T is invalid."))
#     n = [0,0,0,0]               # counter for 00,01,10,11
#     shift = N-b                 # shift desired bit b
#     shift2 = N-b-1              # shift
#     mask = one(T) << shift      # mask everything except b
#     # get the bit b of the first element in A
#     a1 = (A[1] & mask) >> shift
#     for a in A[2:end]                      # mask everything but b and shift
#
#         n += (a & mask) >> shift    # to have either 0x00 or 0x01
#         n[a1 & a2] += 1
#     end
#     return n
# end

bitcount(A::Array{T},b::Int) where {T<:Union{Signed,AbstractFloat}} =
    bitcount(reinterpret.(whichUInt(T),A),b)

function bitcount(A::Array{T}) where {T<:Union{Unsigned,Signed,AbstractFloat}}
    N = 8*sizeof(T)
    n = fill(0,N)
    for b in 1:N
        n[b] = bitcount(A,b)
    end
    return n
end

function bitcountentropy(A::AbstractArray)
    N = prod(size(A))
    p = bitcount(A) / N
    e = [abs(entropy([pi,1-pi],2)) for pi in p]
    return e
end
