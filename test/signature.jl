@testset "SurvivalSignature" begin

    A = zeros(6, 6)
    A[1, [2, 3]] .= 1.0
    A[2, [4, 5]] .= 1.0
    A[3, [4, 6]] .= 1.0
    A[4, [5, 6]] .= 1.0

    types = Dict(1 => [1, 2, 5], 2 => [3, 4, 6])

    φ = s_t_connectivity([1:6;], [1], [5, 6])

    @testset "exactentry" begin
        @test SurvivalSignature.exactentry(CartesianIndex(1, 1), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(1, 2), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(1, 3), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(1, 4), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(2, 1), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(2, 2), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(2, 3), A, types, φ) == 1.0 / 9.0
        @test SurvivalSignature.exactentry(CartesianIndex(2, 4), A, types, φ) == 3.0 / 9.0
        @test SurvivalSignature.exactentry(CartesianIndex(3, 1), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(3, 2), A, types, φ) == 0.0
        @test SurvivalSignature.exactentry(CartesianIndex(3, 3), A, types, φ) == 4.0 / 9.0
        @test SurvivalSignature.exactentry(CartesianIndex(3, 4), A, types, φ) == 6.0 / 9.0
        @test SurvivalSignature.exactentry(CartesianIndex(4, 1), A, types, φ) == 1.0
        @test SurvivalSignature.exactentry(CartesianIndex(4, 2), A, types, φ) == 1.0
        @test SurvivalSignature.exactentry(CartesianIndex(4, 3), A, types, φ) == 1.0
        @test SurvivalSignature.exactentry(CartesianIndex(4, 4), A, types, φ) == 1.0
    end
end
