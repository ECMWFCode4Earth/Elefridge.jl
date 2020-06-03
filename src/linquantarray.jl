using BitIntegers
BitIntegers.@define_integers 24

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
    A::AbstractArray{UInt24,N}
    min::Float64
    max::Float64
end

LinQuantArray = Union{  LinQuant8Array,
                        LinQuant16Array,
                        LinQuant24Array}

Base.size(QA::LinQuantArray) = size(QA.A)
Base.getindex(QA::LinQuantArray,i...) = getindex(QA.A,i...)
# Base.setindex!(QA::ScaledArray,x::Real,i...) = setindex!(SA.A,x,i...)

function whichUInt(n::Integer)
    n == 8 && return UInt8
    n == 16 && return UInt16
    n == 24 && return UInt24
    throw(error("Only n=8,16,24 supported."))
end

function LinQuantization(n::Integer,A::AbstractArray)
    Amin = Float64(minimum(A))
    Amax = Float64(maximum(A))

    Δ = (Amax-Amin)/(2^n-1)
    bounds = Array{Float64,1}(undef,2^n+1)
    bounds[1] = Amin
    bounds[2] = Amin + Δ/2
    bounds[end] = Amax
    for i in 2:2^n
        bounds[i] = bounds[i-1]+Δ
    end

    s = size(A)
    T = whichUInt(n)
    Q = Array{T,length(s)}(undef,s...)

    for i in eachindex(Q)
        Q[i] = findFirstSmaller(Float64(A[i]),bounds)-1
    end

    return Q,Amin,Amax
end

LinQuant8Array(A::AbstractArray) = LinQuant8Array(LinQuantization(8,A)...)
LinQuant16Array(A::AbstractArray) = LinQuant16Array(LinQuantization(16,A)...)
LinQuant24Array(A::AbstractArray) = LinQuant24Array(LinQuantization(24,A)...)

function Array{T}(n::Integer,Q::LinQuantArray) where T
    s = size(Q)
    A = Array{T,length(s)}(undef,s...)
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

Array(Q::LinQuant8Array) = Array{Float32}(8,Q)
Array(Q::LinQuant16Array) = Array{Float32}(16,Q)
Array(Q::LinQuant24Array) = Array{Float32}(24,Q)
