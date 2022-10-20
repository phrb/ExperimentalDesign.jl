@testset "paley" begin
    M = paley(Matrix{Int}(undef, 8, 8))

    Mt = [-1   1   1  -1   1   1   1  -1
 1   1  -1   1   1   1  -1  -1
 1  -1   1   1   1  -1  -1   1
-1   1   1   1  -1  -1   1   1
 1   1   1  -1  -1   1   1  -1
 1   1  -1  -1   1   1  -1   1
 1  -1  -1   1   1  -1   1   1
-1  -1   1   1  -1   1   1   1]

    @test M == Mt
end
