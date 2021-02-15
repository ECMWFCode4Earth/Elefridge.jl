"""Calculte the Continuous Ranked Probability Score (CRPS) of an ensemble array wrt
a true observation o. In-place version that will sort the array."""
function CRPS!(x::Array{T,1},o::T) where T
    sort!(x)
    n = length(x)
    dp = 1/n                                # probability increment per xi in x

    # find index m to split integration over x into heaviside=0 or 1
    m = findfirst(xi -> xi >= o,x)   
    
    # integration from o to x if o is outside of the range of x
    s = m==1 ? x[1]-o : 0.0                 

    # integrate over the heaviside=0 part
    for i in 1:(m == nothing ? n-1 : m)
        s += (i*dp)^2*(x[i+1]-x[i])
    end

    # integrate over the heaviside=0 part
    for i in (m == nothing ? n : m+1):n-1
        s += (i*dp-1)^2*(x[i+1]-x[i])
    end

    s += (m == nothing) ? o-x[end] : 0.0

    return s
end

"""Calculte the Continuous Ranked Probability Score (CRPS) of an ensemble array wrt
a true observation o."""
CRPS(x::Array{T,1},o::T) where T = CRPS!(copy(x),o)
CRPS(x::Array{T,1},o::Real) where T = CRPS!(copy(x),T(o))

