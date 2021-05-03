using Test
using NanoporeMultiEventalign

@testset "dtw.jl" begin
    @testset "nanopore_dtw" begin
        @test nanopore_dtw(
        convert(Vector{Float32}, [1,7,1,6,0,1,6,1]),
        convert(Vector{Float32}, [1,7,1,6,0,1,6,1]),
        "../models/r9.4_70bps.u_to_t_rna.5mer.template.model"
        )[1] == 0.0

        @test nanopore_dtw(
        convert(Vector{Float32}, [1,7,1,6,0,1,6,1]),
        convert(Vector{Float32}, [1,6,1,0,6,1,7,1]),
        "../models/r9.4_70bps.u_to_t_rna.5mer.template.model"
        )[1] != 0.0

        @test nanopore_dtw( # tests the chengepoints detection
        convert(Vector{Float32}, [1,2,3,4,5,6,7]),
        convert(Vector{Float32}, [7,6,5,4,3,2,1]),
        "../models/r9.4_70bps.u_to_t_rna.5mer.template.model"
        )[1] == 0.0

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
    @testset "bhattacharyya" begin
        @test bhattacharyya(kmerdist(1,1),kmerdist(1,1)) == 0.0
        @test bhattacharyya(kmerdist(1,2),kmerdist(3,2)) == 0.029445758914095864
    end
end
