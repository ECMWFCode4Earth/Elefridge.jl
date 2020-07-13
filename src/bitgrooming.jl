"""Creates a UInt32-mask for the trailing non-significant bits of a
Float32 number. `nsb` are the number of significant bits in the mantissa.
E.g. mask(3) returns `00000000000011111111111111111111`,
such that all but the first 3 significant bits can be masked."""
mask(nsb::Integer) = UInt32(2^(23-nsb)-1)

"""Shave trailing bits of a Float32 number to zero.
Mask is UInt32 with 1 for the shaved bits, 0 for the retained bits."""
function shave(x::Float32,mask::UInt32)
    ui = reinterpret(UInt32,x)
    ui &= mask
    return reinterpret(Float32,ui)
end

"""Shave trailing bits of a Float32 number to zero.
Providing `nsb` the number of retained significant bits, a mask is created
and applied."""
shave(x::Float32,nsb::Integer) = shave(x,~mask(nsb))

"""Shave trailing bits of a Float32 number to zero.
In case no `sb` argument is applied for `shave`, shave 16 bits, retain 7."""
shave(x::Float32) = shave(x,7)

"""Shave trailing bits of a Float32 array to zero.
Creates the shave-mask only once and applies it to every element in `X`."""
shave(X::AbstractArray{Float32},nsb::Integer) = shave.(X,~mask(nsb))

"""Set trailing bits of a Float32 number to one.
Provided a UInt32 mask with 1 for bits to be set to one, and 0 else."""
function set_one(x::Float32,mask::UInt32)
    ui = reinterpret(UInt32,x)
    ui |= mask
    return reinterpret(Float32,ui)
end

"""Set trailing bits of Float32 number to one, given `nsb` number of significant
bits retained. A mask is created and applied."""
set_one(x::Float32,nsb::Integer) = set_one(x,mask(nsb))

"""Set trailing bits of a Float32 number to one.
In case no `sb` argument is applied for `set_one`, set 16 bits, retain 7."""
set_one(x::Float32) = set_one(x,7)

"""Set trailing bits of a Float32 number to one.
Creates the setting-mask only once and applies it to every element in `X`."""
set_one(X::AbstractArray{Float32},nsb::Integer) = set_one.(X,mask(nsb))

"""Bit-grooming. Alternatingly apply bit-shaving and setting to a Float32 array."""
function groom(X::AbstractArray{Float32},nsb::Integer)
    Y = similar(X)      # preallocate output of same size and type
    mask1 = mask(nsb)   # mask for setting
    mask0 = ~mask1      # mask for shaving
    @inbounds for i in 1:2:length(X)-1
        Y[i] = shave(X[i],mask0)            # every second element is shaved
        Y[i+1] = set_one(X[i+1],mask1)      # every other 2nd element is set
    end
    return Y
end

"""Calculates an integeer as argument for a bitshift operation
required to move the least significant bit (after rounding) to
the last bit."""
shift(nsb::Integer) = 23-nsb

"""Creates a mask for bit-setting given `nsb` bits to be retained in the
significand."""
setmask(nsb::Integer) = 0x003f_ffff >> nsb

"""Round to nearest for floating-point arithmetic, using only integer
arithmetic. `setmask`,`shift`,`shavemask` have to be provided that depend
on the number of significant bits that will be retained."""
function Base.round(x::Float32,
                    setmask::UInt32,
                    shift::Integer,
                    shavemask::UInt32)
    ui = reinterpret(UInt32,x)
    ui += setmask + ((ui >> shift) & 0x0000_0001)
    return reinterpret(Float32,ui & shavemask)
end

"""Round to nearest for Float32, given `nsb` number of signifcant bits, that
will be retained. E.g. round(x,7) will round the trailing 16 bits and retain
the 7 significant bits (which might be subject to change by a carry bit)."""
Base.round(x::Float32,nsb::Integer) = round(x,setmask(nsb),shift(nsb),~mask(nsb))

"""Round to nearest for a Float32 array `X`. The bit-masks are only created once
and then applied to every element in `X`."""
function Base.round(X::AbstractArray{Float32},nsb::Integer)
    semask = setmask(nsb)
    s = shift(nsb)
    shmask = ~mask(nsb)
    return round.(X,semask,s,shmask)
end

"""Number of significant bits `nsb` given the number of significant digits `nsd`."""
nsb(nsd::Integer) = Integer(ceil(log(10)/log(2)*nsd))
