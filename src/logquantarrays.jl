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
    A::AbstractArray{UInt32,N}
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

    logmin = log10(Float64(minimum(A)))
    logmax = log10(Float64(maximum(A)))

    Δ = (logmax-logmin)/(2^n-1)
    bounds = Array{Float64,1}(undef,2^n+1)
    bounds[1] = logmin          # bounds in log space
    bounds[2] = logmin + Δ/2
    bounds[end] = logmax
    for i in 2:2^n
        bounds[i] = bounds[i-1]+Δ
    end
    bounds = 10.0 .^ bounds     # only now convert to lin space

    s = size(A)
    T = whichUInt(n)
    Q = Array{T,length(s)}(undef,s...)

    for i in eachindex(Q)
        Q[i] = findFirstSmaller(Float64(A[i]),bounds)-1
    end

    return Q,logmin,logmax
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
        A[i] = ten^(Qlogmin + Q[i]*Δ)
    end
    return A
end

Array{T}(Q::LogQuant8Array) where T = Array{T}(8,Q)
Array{T}(Q::LogQuant16Array) where T = Array{T}(16,Q)
Array{T}(Q::LogQuant24Array) where T = Array{T}(24,Q)

Array(Q::LogQuant8Array) = Array{Float32}(8,Q)
Array(Q::LogQuant16Array) = Array{Float32}(16,Q)
Array(Q::LogQuant24Array) = Array{Float32}(24,Q)
