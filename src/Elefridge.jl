module Elefridge

    export LinQuant8Array, LinQuant16Array, LinQuant24Array,
        LogQuant8Array, LogQuant16Array, LogQuant24Array,
        bitentropy, shave, set_one, groom


    using StatsBase
    using BitIntegers
    BitIntegers.@define_integers 24

    include("linquantarrays.jl")
    include("logquantarrays.jl")
    include("bitentropy.jl")
    include("bitgrooming.jl")

end
