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

whichUInt(::Type{T}) where T = whichUInt(sizeof(T)*8)
