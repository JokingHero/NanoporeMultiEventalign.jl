using NanoporeMultiEventalign
using Test

@testset "dtw.jl" begin
    @testset "dtw" begin
        @test dtw(true)
        @test !dtw(false)
    end
end

@testset "persistence.jl" begin
    @testset "loadnanoporefast5" begin
        @test loadnanoporefast5("test/sample_data/fast5/control/0/0a4c39b4-14c0-4d49-bce9-81486627062c.fast5")[6] == 76.81596
    end
end
