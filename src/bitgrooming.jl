mask(nsb::Integer) = UInt32(2^(23-nsb)-1)

function shave(x::Float32,mask::UInt32)
    ui = reinterpret(UInt32,x)
    ui &= mask
    return reinterpret(Float32,ui)
end

shave(x::Float32,nsb::Integer) = shave(x,~mask(nsb))
shave(x::Float32) = shave(x,7)
shave(X::AbstractArray{Float32},nsb::Integer) = shave.(X,~mask(nsb))

function set_one(x::Float32,mask::UInt32)
    ui = reinterpret(UInt32,x)
    ui |= mask
    return reinterpret(Float32,ui)
end

set_one(x::Float32,nsb::Integer) = set_one(x,mask(nsb))
set_one(x::Float32) = set_one(x,7)
set_one(X::AbstractArray{Float32},nsb::Integer) = set_one.(X,mask(nsb))

function groom(X::AbstractArray{Float32},nsb::Integer)
    Y = similar(X)
    mask1 = mask(nsb)
    mask0 = ~mask1
    @inbounds for i in 1:2:length(X)-1
        Y[i] = shave(X[i],mask0)
        Y[i+1] = set_one(X[i+1],mask1)
    end
    return Y
end

shift(nsb::Integer) = 23-nsb
setmask(nsb::Integer) = 0x003f_ffff >> nsb

function Base.round(x::Float32,
                    setmask::UInt32,
                    shift::Integer,
                    shavemask::UInt32)
    ui = reinterpret(UInt32,x)
    ui += setmask + ((ui >> shift) & 0x0000_0001)
    return reinterpret(Float32,ui & shavemask)
end

Base.round(x::Float32,nsb::Integer) = round(x,setmask(nsb),shift(nsb),~mask(nsb))

function Base.round(X::AbstractArray{Float32},nsb::Integer)
    semask = setmask(nsb)
    s = shift(nsb)
    shmask = ~mask(nsb)
    return round.(X,semask,s,shmask)
end

nsb(nsd::Integer) = Integer(ceil(log(10)/log(2)*nsd))
