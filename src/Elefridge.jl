module Elefridge

    export CRPS, CRPS!

    using StatsBase

    include("helpers.jl")
    include("nzextrema.jl")
    include("linquantarrays.jl")
    include("logquantarrays.jl")
    include("bitentropy.jl")
    include("bitgrooming.jl")
    include("bitrounding.jl")
    include("bitstring.jl")
    include("bitinformation.jl")
    include("bittranspose.jl")
    include("xor_delta.jl")
    include("crps.jl")
end
