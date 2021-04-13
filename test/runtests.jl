using Test
using Elefridge

@testset "CRPS tests" begin

    # Gaussian around true 0
    c = 0.2336949772551091
    @test isapprox(CRPS(randn(100_000),0.0),c,rtol=1e-1)

    # for true = 1,-1 the score needs to increase
    @test CRPS(randn(100_000),1) > c
    @test CRPS(randn(100_000),-1) > c
end
