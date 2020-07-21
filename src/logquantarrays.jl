struct LogQuantArray{T,N} <: AbstractArray{Unsigned, N}
    A::Array{T,N}
    min::Float64
    max::Float64
end

Base.size(QA::LogQuantArray) = size(QA.A)
Base.getindex(QA::LogQuantArray,i...) = getindex(QA.A,i...)
Base.eltype(Q::LogQuantArray{T,N}) where {T,N} = T

function LogQuantization(n::Integer,A::Array{T2,N}) where {T2,N}

    any(A .<= zero(eltype(A))) && throw(DomainError(
                    "LogQuantization only for positive arguments."))

    logmin = log(minimum(A))
    logmax = log(maximum(A))
    Δ = (2^n-1)/(logmax-logmin)

    T = whichUInt(n)
    Q = similar(A,T)

    @. @views Q = T(round((log(A)-logmin)*Δ))

    return LogQuantArray{T,N}(Q,Float64(logmin),Float64(logmax))
end

LogQuant8Array(A::Array{T,N}) where {T,N} = LogQuantization(8,A)
LogQuant16Array(A::Array{T,N}) where {T,N} = LogQuantization(16,A)
LogQuant24Array(A::Array{T,N}) where {T,N} = LogQuantization(24,A)

function Array{T}(n::Integer,Q::LogQuantArray) where T
    Qlogmin = T(Q.min)
    Qlogmax = T(Q.max)
    Δ = (Qlogmax-Qlogmin)/(2^n-1)

    A = similar(Q,T)

    @inbounds for i in eachindex(A)
        A[i] = exp(Qlogmin + Q[i]*Δ)
    end

    return A
end

Array{T}(Q::LogQuantArray{UInt8,N}) where {T,N} = Array{T}(8,Q)
Array{T}(Q::LogQuantArray{UInt16,N}) where {T,N} = Array{T}(16,Q)
Array{T}(Q::LogQuantArray{UInt24,N}) where {T,N} = Array{T}(24,Q)

Array(Q::LogQuantArray{UInt8,N}) where N = Array{Float32}(8,Q)
Array(Q::LogQuantArray{UInt16,N}) where N = Array{Float32}(16,Q)
Array(Q::LogQuantArray{UInt24,N}) where N = Array{Float32}(24,Q)
