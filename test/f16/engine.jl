using Test
using CSV
using DataFrames
using FlightMechanicsSimulator
using FlightMechanicsUtils


@testset "engine.jl" begin

    df = DataFrame(CSV.File("data/tgear.csv"))
    for case in eachrow(df)
        rv1 = F16Stevens.tgear(case.thtl)
        @test isapprox(rv1, case.tgear, atol = 1.0e-15)
    end

    df = DataFrame(CSV.File("data/pdot.csv"))
    for case in eachrow(df)
        rv1 = F16Stevens.pdot(case.p3, case.p1)
        @test isapprox(rv1, case.pdot, atol = 1.0e-15)
    end

    df = DataFrame(CSV.File("data/rtau.csv"))
    for case in eachrow(df)
        rv1 = F16Stevens.rtau(case.dp)
        @test isapprox(rv1, case.rtau, atol = 1.0e-15)
    end

    df = DataFrame(CSV.File("data/thrust.csv"))
    for case in eachrow(df)
        rv1 = F16Stevens.thrust(case.pow, case.alt * FT2M, case.rmach)
        @test isapprox(rv1 / (LB2KG * F16Stevens.GD * FT2M), case.thrust, atol = 1.0e-10)
    end

end
