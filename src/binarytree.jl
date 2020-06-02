function ispower2(x::Integer)
    while true
        if x == 2
            return true
        elseif x % 2 == 1
            return false
        else
            x = x ÷ 2
        end
    end
    return false
end

function findFirstSmaller(x::Float64,v::Array{Float64,1})

    l = length(v)-1
    @boundscheck ispower2(l) || throw(BoundsError())
    n = Int(log2(l))-1
    idx = l ÷ 2

    if x > v[end-1]     # round to max
        return l
    end

    # binary tree search
    for i in 1:n
        if x < v[idx]
            idx -= 2^(n-i)
        else
            idx += 2^(n-i)
        end
    end

    # split off the i = n+1 case
    if x >= v[idx]
        idx += 1
    end

    return idx-1
end
