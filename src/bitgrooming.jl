function shave(x::Float32,mask::UInt32)
    ui = reinterpret(UInt32,x)
    ui &= mask
    return reinterpret(Float32,ui)
end

shave(x::Float32,nsb::Integer) = shave(x,~UInt32(2^(23-nsb)-1))

function shave(X::AbstractArray{Float32},nsb::Integer)
    mask = ~UInt32(2^(23-nsb)-1)
    return shave.(X,mask)
end

function set_one(x::Float32,mask::UInt32)
    ui = reinterpret(UInt32,x)
    ui |= mask
    return reinterpret(Float32,ui)
end

set_one(x::Float32,nsb::Integer) = set_one(x,UInt32(2^(23-nsb)-1))

function set_one(X::AbstractArray{Float32},nsb::Integer)
    mask = UInt32(2^(23-nsb)-1)
    return set_one.(X,mask)
end

function groom(X::AbstractArray{Float32},nsb::Integer)
    Y = similar(X)
    mask = UInt32(2^(23-nsb)-1)
    mask_inv = ~mask
    @inbounds for i in 1:2:length(X)-1
        Y[i] = shave(X[i],mask_inv)
        Y[i+1] = set_one(X[i],mask)
    end
    return Y
end

function Base.round(x::Float32,nsb::Integer)
    ui = reinterpret(UInt32,x)
    ui += 0x0000_7fff + ((ui >> 16) & 0x0000_0001)
    return reinterpret(Float32,ui & 0xffff0000)
end

Base.round(X::AbstractArray{Float32},nsb::Integer) = round.(X,nsb)

nsb(nsd::Integer) = Integer(ceil(log(10)/log(2)*nsd))

# Some tests
using Statistics

A = Array{Float32}(1:1e-6:2-1e-6)
Ao = set_one.(A,8)
As = shave.(A,8)
Ar = round.(A,sigdigits=9,base=2)
Ag = groom(A,8)

# check bits that are always zero
zero_check = UInt32(0)
for x in Ar
    global zero_check |= reinterpret(UInt32,x)
end
bitstring(zero_check)

# bias?
mean(As)-mean(A)
mean(Ao)-mean(A)
mean(Ar)-mean(A)
mean(Ag)-mean(A)

# mean error
mean(As-A)
mean(Ao-A)
mean(Ar-A)
mean(Ag-A)

# mean absolute error
mean(abs.(As-A))
mean(abs.(Ao-A))
mean(abs.(Ar-A))
mean(abs.(Ag-A))

# max absolute error
maximum(abs.(As-A))
maximum(abs.(Ao-A))
maximum(abs.(Ar-A))
maximum(abs.(Ag-A))

# max decimal error
maximum(abs.(log10.(As./A)))
maximum(abs.(log10.(Ao./A)))
maximum(abs.(log10.(Ar./A)))
maximum(abs.(log10.(Ag./A)))
