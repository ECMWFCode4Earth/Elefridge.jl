"""Return minimum, ignoring negatives and zeroes."""
function minpos(A::AbstractArray{T}) where T
    o = zero(T)
    foldl((x,y) -> y > o ? min(x,y) : x, A; init=typemax(T))
end
