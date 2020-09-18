@testset "XOR reversibility UInt" begin
    for T in (UInt8,UInt16,UInt32,UInt64)
        A = rand(T,1000)
        @test A == unxor_delta(xor_delta(A))
        @test A == xor_delta(unxor_delta(A))

        B = copy(A)
        xor_delta!(A)
        unxor_delta!(A)
        @test B == A

        unxor_delta!(A)
        xor_delta!(A)
        @test B == A
    end
end

@testset "XOR reversibility Float" begin
    for T in (Float32,Float64)
        A = rand(T,1000)
        @test A == unxor_delta(xor_delta(A))
        @test A == xor_delta(unxor_delta(A))
    end
end
