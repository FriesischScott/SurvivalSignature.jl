@testset "Percolation" begin

    @test isapprox(percolation(gridnetwork(5, 5)), 0.574468, rtol = 1e-5)
    @test isapprox(percolation(gridnetwork(5, 6)), 0.584746, rtol = 1e-5)
    @test isapprox(percolation(gridnetwork(6, 6)), 0.594595, rtol = 1e-5)
end
