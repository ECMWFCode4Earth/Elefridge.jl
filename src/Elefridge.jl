module Elefridge

    export  LinQuantArray, LogQuantArray, 
        LinQuant8Array, LinQuant16Array, LinQuant24Array,
        LogQuant8Array, LogQuant16Array, LogQuant24Array,
        bitentropy, shave, set_one, groom


    using StatsBase
    using BitIntegers
    BitIntegers.@define_integers 24

    export UInt24

    include("linquantarrays.jl")
    include("logquantarrays.jl")
    include("bitentropy.jl")
    include("bitgrooming.jl")
    include("bitrounding.jl")

end
