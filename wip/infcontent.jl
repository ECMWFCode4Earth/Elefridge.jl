function bitcount(A::Array{UInt8,1})
    n = fill(0,8)
    mask = 0x01
    for a in A
        for i in 0:7
            n[i+1] += (a & (mask << i)) >> i
        end
    end
    return n
end

function bitcounts(a::Array{UInt8,1})
    n = fill(0,8)
    for b in 1:8
        n[b] = bitcount(a,b)
    end
    n
end

function bitcount(A::Array{UInt8,1},b::Int)
    n = 0
    shift = b-1
    mask = 0x01 << shift
    for a in A
        n += (a & mask) >> shift
    end
    return n
end
