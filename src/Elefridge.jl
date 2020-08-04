module Elefridge

    export  LinQuantArray, LogQuantArray,
        LinQuant8Array, LinQuant16Array, LinQuant24Array,
        LogQuant8Array, LogQuant16Array, LogQuant24Array,
        bitentropy, shave, set_one, groom

    export bitstring


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

end
