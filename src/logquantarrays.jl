struct LogQuantArray{T,N} <: AbstractArray{Unsigned, N}
    A::Array{T,N}
    min::Float64
    max::Float64
end

Base.size(QA::LogQuantArray) = size(QA.A)
Base.getindex(QA::LogQuantArray,i...) = getindex(QA.A,i...)
Base.eltype(Q::LogQuantArray{T,N}) where {T,N} = T

function LogQuantization(   ::Type{T},
                            A::AbstractArray,
                            round_nearest_in::Symbol=:linspace) where T

    (any(A .< zero(eltype(A))) || ~all(isfinite.(A))) &&
        throw(DomainError("LogQuantization only for positive&zero entries."))

    # min/max of non-zero entries
    mi = Float64(minpos(A))
    logmin = log(Float64(mi))
    logmax = log(Float64(maximum(A)))

    # throw error in case the range is zero.
    logmin == logmax && throw(DomainError("Data range is zero."))

    # range of values in log-space
    # map min to 1 and max to ff..., reserve 0 for 0.
    Δ = (2^(sizeof(T)*8)-2)/(logmax-logmin)

    # shift to round-to-nearest in lin or log-space
    if round_nearest_in == :linspace
        c = 1/2 - Δ*log(mi*(exp(1/Δ)+1)/2)
    elseif round_nearest_in == :logspace
        c = -logmin*Δ
    else
        throw(ArgumentError("Round-to-nearest either :linspace or :logspace"))
    end

    # preallocate output
    Q = similar(A,T)

    for i in eachindex(A)
        if iszero(A[i]) # store as 0x00...
            Q[i] = zero(T)
        else    # positive numbers: convert to logpacking in 0x1-0xff..
            Q[i] = T(round(c + Δ*log(Float64(A[i]))))+one(T)
        end
    end

    return LogQuantArray{T,ndims(Q)}(Q,Float64(logmin),Float64(logmax))
end

function LogQuant8Array(A::Array{T,N},rn::Symbol=:linspace) where {T,N}
    LogQuantization(whichUInt(8),A,rn)
end

function LogQuant16Array(A::Array{T,N},rn::Symbol=:linspace) where {T,N}
     LogQuantization(whichUInt(16),A,rn)
end

function LogQuant24Array(A::Array{T,N},rn::Symbol=:linspace) where {T,N}
     LogQuantization(whichUInt(24),A,rn)
end

function Array{T}(n::Integer,Q::LogQuantArray) where T
    Qlogmin = T(Q.min)
    Qlogmax = T(Q.max)
    Δ = (Qlogmax-Qlogmin)/(2^n-2) # -2 as 0x00.. is reserved for 0

    A = similar(Q,T)

    @inbounds for i in eachindex(A)
        if iszero(Q[i])
            A[i] = zero(T)
        else
            A[i] = exp(Qlogmin + (Q[i]-1)*Δ)
        end
    end

    return A
end

Array{T}(Q::LogQuantArray{UInt8,N}) where {T,N} = Array{T}(8,Q)
Array{T}(Q::LogQuantArray{UInt16,N}) where {T,N} = Array{T}(16,Q)
Array{T}(Q::LogQuantArray{UInt24,N}) where {T,N} = Array{T}(24,Q)

Array(Q::LogQuantArray{UInt8,N}) where N = Array{Float32}(8,Q)
Array(Q::LogQuantArray{UInt16,N}) where N = Array{Float32}(16,Q)
Array(Q::LogQuantArray{UInt24,N}) where N = Array{Float32}(24,Q)
