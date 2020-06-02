struct LinQuant8Array{N} <: AbstractArray{UInt8, N}
    A::AbstractArray{UInt8,N}
    min::Float64
    max::Float64
end

struct LinQuant16Array{N} <: AbstractArray{UInt16, N}
    A::AbstractArray{UInt16,N}
    min::Float64
    max::Float64
end

struct LinQuant24Array{N} <: AbstractArray{UInt32, N}
    A::AbstractArray{UInt32,N}
    min::Float64
    max::Float64
end

struct LinQuant32Array{N} <: AbstractArray{UInt32, N}
    A::AbstractArray{UInt32,N}
    min::Float64
    max::Float64
end

LinQuantArray = Union{  LinQuant8Array,
                        LinQuant16Array,
                        LinQuant24Array,
                        LinQuant32Array}

Base.size(QA::LinQuantArray) = size(QA.A)
Base.getindex(QA::LinQuantArray,i...) = getindex(QA.A,i...)
# Base.setindex!(QA::ScaledArray,x::Real,i...) = setindex!(SA.A,x,i...)

function whichUInt(n::Integer)
    n == 8 && return UInt8
    n == 16 && return UInt16
    n == 24 && return UInt32
    n == 32 && return UInt32
    throw(error("Only n=6,16,24,32 supported."))
end

function LinQuantArray(n::Integer,A::AbstractArray)
    Amin = Float64(minimum(A))
    Amax = Float64(maximum(A))

    Δ = (Amax-Amin)/(2^n-1)
    quants = Array(Amin:Δ:Amax)
    bounds = vcat(Amin,(quants[1:end-1]+quants[2:end])/2,Amax)

    s = size(A)
    T = whichUInt(n)
    Q = Array{T,length(s}(undef,s...)

    for i in eachindex(Q)
        Q[i] = findFirstSmaller(Float64(A[i]),bounds)-1
    end

    return Q,Amin,Amax
end

LinQuant8Array(A::AbstractArray) = LinQuant8Array(LinQuantArray(8,A)...)
LinQuant16Array(A::AbstractArray) = LinQuant16Array(LinQuantArray(16,A)...)
LinQuant24Array(A::AbstractArray) = LinQuant24Array(LinQuantArray(24,A)...)
LinQuant32Array(A::AbstractArray) = LinQuant32Array(LinQuantArray(32,A)...)

function Array{T}(n::Integer,Q::LinQuantArray{N}) where {T,N}
    A = Array{T,N}(undef,size(Q)...)
    Qmin = T(Q.min)
    Qmax = T(Q.max)
    Δ = (Qmax-Qmin)/(2^n-1)

    for i in eachindex(A)
        A[i] = Qmin + Q[i]*Δ
    end
    return A
end

Array{T}(Q::LinQuant8Array) where T = Array{T}(8,Q)
Array{T}(Q::LinQuant16Array) where T = Array{T}(16,Q)
Array{T}(Q::LinQuant24Array) where T = Array{T}(24,Q)
Array{T}(Q::LinQuant32Array) where T = Array{T}(32,Q)
