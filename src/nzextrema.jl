function nzminimum(A::AbstractArray)
    o = zero(eltype(A))
    foldl((x,y) -> y > o ? min(x,y) : x, A; init=Inf)
end
