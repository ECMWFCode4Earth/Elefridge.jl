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
