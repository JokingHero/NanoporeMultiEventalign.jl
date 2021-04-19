using Test
using NanoporeMultiEventalign

@testset "dtw.jl" begin
    @testset "dtw" begin
        #@test dtw(true)
        #@test !dtw(false)
    end
end

@testset "persistence.jl" begin
    @testset "loadnanoporefast5" begin
        @test loadnanoporefast5(
            "sample_data/fast5/control/0/" *
            "0a4c39b4-14c0-4d49-bce9-81486627062c.fast5")[6] == Float32(76.81596)
    end
    @testset "loadfasta" begin
        @test convert(String, loadfasta(
            "sample_data/reference.fa")[1]) == ""*
            "GAATACAAGCTACTTGTTCTTTTTGCAGGATCCCATCGATTCGAATTCAAGGCCTCTCGAGCCTCTA"*
            "GAACTATAGTGAGTCGTATTACGTAGATCCAGACATGATAAGATACATTGATGAGTTTGGACAAACC"*
            "ACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAA"*
            "CCATTATAAGCTGCAATAAACAAGTTAACAACAACAATTGCATTCATTTTATGTTTCAGGTTCAGGG"*
            "GGAGGTGTGGGAGGTTTTTTAATTCGCAAAAAAAAAAGGGAGAGCGACGAGAAGAAGAGGCTTGTAT"*
            "GGACGATGTCTTCCACAACAGAAGACATCGTCCATACAAGCAAAAAA"
    end
end

@testset "utils.jl" begin
    @testset "Bhattacharyya" begin
        @test bhattacharyya(kmerdist(1,1),kmerdist(1,1)) == 0.0
        @test bhattacharyya(kmerdist(1,2),kmerdist(3,2)) == 0.029445758914095864
    end
end
