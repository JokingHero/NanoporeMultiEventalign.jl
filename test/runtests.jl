using Test
using NanoporeMultiEventalign

@testset "dtw.jl" begin
    @testset "dtw" begin
        @test dtw(true)
        @test !dtw(false)
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
