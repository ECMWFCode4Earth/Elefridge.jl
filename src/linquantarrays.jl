struct LinQuantArray{T,N} <: AbstractArray{Unsigned,N}
    A::Array{T,N}
    min::Float64
    max::Float64
end

Base.size(QA::LinQuantArray) = size(QA.A)
Base.getindex(QA::LinQuantArray,i...) = getindex(QA.A,i...)

function whichUInt(n::Integer)
    n == 8 && return UInt8
    n == 16 && return UInt16
    n == 24 && return UInt24
    throw(error("Only n=8,16,24 supported."))
end

function LinQuantization(n::Integer,A::Array{T2,N}) where {T2,N}
    Amin = Float64(minimum(A)) # some precision issue when using Float32
    Amax = Float64(maximum(A))
    Δ = (2^n-1)/(Amax-Amin)

    T = whichUInt(n)
    Q = similar(A,T)

    @. @views Q = T(round((A-Amin)*Δ))

    return LinQuantArray{T,N}(Q,Amin,Amax)
end

LinQuant8Array(A::Array{T,N}) where {T,N} = LinQuantization(8,A)
LinQuant16Array(A::Array{T,N}) where {T,N} = LinQuantization(16,A)
LinQuant24Array(A::Array{T,N}) where {T,N} = LinQuantization(24,A)

function Array{T}(n::Integer,Q::LinQuantArray) where T
    Qmin = T(Q.min)
    Qmax = T(Q.max)
    Δ = (Qmax-Qmin)/(2^n-1)

    A = similar(Q,T)

    @inbounds for i in eachindex(A)
        A[i] = Qmin + Q[i]*Δ
    end

    return A
end

Array{T}(Q::LinQuantArray{UInt8,N}) where {T,N} = Array{T}(8,Q)
Array{T}(Q::LinQuantArray{UInt16,N}) where {T,N} = Array{T}(16,Q)
Array{T}(Q::LinQuantArray{UInt24,N}) where {T,N} = Array{T}(24,Q)

Array(Q::LinQuantArray{UInt8,N}) where N = Array{Float32}(8,Q)
Array(Q::LinQuantArray{UInt16,N}) where N = Array{Float32}(16,Q)
Array(Q::LinQuantArray{UInt24,N}) where N = Array{Float32}(24,Q)
