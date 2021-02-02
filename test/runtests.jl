using NanoporeMultiEventalign
using Test

@testset "dtw.jl" begin
    @testset "dtw" begin
        @test dtw(true)
        @test !dtw(false)
    end
end
