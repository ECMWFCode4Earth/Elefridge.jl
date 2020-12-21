module Elefridge

    export  LinQuantArray, LogQuantArray,
        LinQuant8Array, LinQuant16Array, LinQuant24Array,
        LogQuant8Array, LogQuant16Array, LogQuant24Array,
        shave, set_one, groom, halfshave, minpos, round!

    export bitstring, bitentropy, bitcount, bitpaircount,
        bitinformation, signed_exponent, signed_exponent!,
        bittranspose, bitbacktranspose, xor_delta,
        unxor_delta


    using StatsBase
    using BitIntegers
    BitIntegers.@define_integers 24
    BitIntegers.@define_integers 40
    BitIntegers.@define_integers 48
    BitIntegers.@define_integers 56

    export UInt24, UInt40, UInt48, UInt56

    include("helpers.jl")
    include("nzextrema.jl")
    include("linquantarrays.jl")
    include("logquantarrays.jl")
    include("bitentropy.jl")
    include("bitgrooming.jl")
    include("bitrounding.jl")
    include("bitstring.jl")
    include("infcontent.jl")
    include("bittranspose.jl")
    include("xor_delta.jl")

end
