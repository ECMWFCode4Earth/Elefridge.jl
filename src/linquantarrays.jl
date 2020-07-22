struct LinQuantArray{T,N} <: AbstractArray{Unsigned,N}
    A::Array{T,N}
    min::Float64
    max::Float64
end

Base.size(QA::LinQuantArray) = size(QA.A)
Base.getindex(QA::LinQuantArray,i...) = getindex(QA.A,i...)
Base.eltype(Q::LinQuantArray{T,N}) where {T,N} = T

function whichUInt(n::Integer)
    n == 8 && return UInt8
    n == 16 && return UInt16
    n == 24 && return UInt24
    n == 32 && return UInt32
    n == 40 && return UInt40
    n == 48 && return UInt48
    n == 56 && return UInt56
    n == 64 && return UInt64
    throw(error("Only n=8,16,24,32,40,48,56,64 supported."))
end

function LinQuantization(n::Integer,A::Array{T2,N}) where {T2,N}

    # determine the range
    Amin = Float64(minimum(A))
    Amax = Float64(maximum(A))
    Δ = (2^n-1)/(Amax-Amin)     # range of values in linear space

    T = whichUInt(n)
    Q = similar(A,T)

    # map minimum to 0, maximum to ff
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
