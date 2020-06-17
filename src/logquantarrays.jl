struct LogQuant8Array{N} <: AbstractArray{UInt8, N}
    A::AbstractArray{UInt8,N}
    min::Float64
    max::Float64
end

struct LogQuant16Array{N} <: AbstractArray{UInt16, N}
    A::AbstractArray{UInt16,N}
    min::Float64
    max::Float64
end

struct LogQuant24Array{N} <: AbstractArray{UInt32, N}
    A::AbstractArray{UInt24,N}
    min::Float64
    max::Float64
end

LogQuantArray = Union{  LogQuant8Array,
                        LogQuant16Array,
                        LogQuant24Array}

Base.size(QA::LogQuantArray) = size(QA.A)
Base.getindex(QA::LogQuantArray,i...) = getindex(QA.A,i...)

function LogQuantization(n::Integer,A::AbstractArray)

    any(A .<= zero(eltype(A))) && throw(DomainError(
                    "LogQuantization only for positive arguments."))

    logmin = log(minimum(A))
    logmax = log(maximum(A))
    Δ = (2^n-1)/(logmax-logmin)

    s = size(A)
    T = whichUInt(n)
    Q = Array{T,length(s)}(undef,s...)

    for i in eachindex(Q)
        Q[i] = T(round((log(A[i])-logmin)*Δ))
    end

    return Q,Float64(logmin),Float64(logmax)
end

LogQuant8Array(A::AbstractArray) = LogQuant8Array(LogQuantization(8,A)...)
LogQuant16Array(A::AbstractArray) = LogQuant16Array(LogQuantization(16,A)...)
LogQuant24Array(A::AbstractArray) = LogQuant24Array(LogQuantization(24,A)...)

function Array{T}(n::Integer,Q::LogQuantArray) where T
    s = size(Q)
    A = Array{T,length(s)}(undef,s...)
    Qlogmin = T(Q.min)
    Qlogmax = T(Q.max)
    Δ = (Qlogmax-Qlogmin)/(2^n-1)

    ten = T(10)
    for i in eachindex(A)
        A[i] = exp(Qlogmin + Q[i]*Δ)
    end
    return A
end

Array{T}(Q::LogQuant8Array) where T = Array{T}(8,Q)
Array{T}(Q::LogQuant16Array) where T = Array{T}(16,Q)
Array{T}(Q::LogQuant24Array) where T = Array{T}(24,Q)

Array(Q::LogQuant8Array) = Array{Float32}(8,Q)
Array(Q::LogQuant16Array) = Array{Float32}(16,Q)
Array(Q::LogQuant24Array) = Array{Float32}(24,Q)
