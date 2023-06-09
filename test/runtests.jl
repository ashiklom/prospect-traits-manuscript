using Test

using ProspectTraits
using CanopyOptics
using Distributions
using LinearAlgebra
using Unitful

@testset "ProspectTraits" begin
    @testset "Spectrum" begin
        @test_warn "No wavelength units given. Assuming nm." Spectrum(400.0:1.0:2500.0, rand(2101))
        @test_nowarn Spectrum((400.0:1.0:2500.0)u"nm", rand(2101))
        spec = Spectrum((400.0:1.0:2500.0)u"nm", rand(2101, 3))
        optis = createLeafOpticalStruct(spec)
    end

    @testset "Fit PROSPECT" begin

        # True values
        N = 1.4
        Ccab = 40.0
        Ccar = 8.0
        Canth = 4.0
        Cbrown = 0.1
        Cw = 0.01
        Cm = 0.01
        Ccbc = 0.008
        Cprot = 0.002
        ρ = 0.999
        σ = 0.007
        
        # Simulate some fake data
        λ_windows = (399.5:1.0:2500.5)*u"nm"
        n = length(λ_windows) - 1
        opti_c = createLeafOpticalStruct(λ_windows; prospect_version="pro")
        true_pro = ProspectTraits.prospectpro(opti_c, N, Ccab, Ccar, Canth, Cbrown,
            Cw, Cprot, Ccbc)
        H = abs.(hcat([collect(1:n) .- i for i in (1:n)]...))
        Ω = ρ.^H
        dist = MvNormal(true_pro, σ^2 * Ω)
        λ = (400.0:1.0:2500.0)*u"nm"
        obs_refl = rand(dist, 3)
        obs = Spectrum(λ, obs_refl)

        @testset "PROSPECT PRO" begin
            r_pro = fit_prospect(obs, 500; version = "pro")
            means = mean(r_pro)
            for sym in (:N, :Ccab, :Ccar, :Canth, :Cbrown, :Cw, :Cprot, :Ccbc,
                    :ρ, :σ)
                @info String(sym)
                @test isapprox(means[sym, :mean], eval(sym), rtol = 0.1)
            end
        end

        @testset "PROSPECT D" begin
            r_d = fit_prospect(obs, 500; version = "d")
            means = mean(r_d)
            for sym in (:N, :Ccab, :Ccar, :Canth, :Cbrown, :Cw, :Cm,
                    :ρ, :σ)
                @info String(sym)
                @test isapprox(means[sym, :mean], eval(sym), rtol = 0.1)
            end
        end

        @testset "PROSPECT 5B" begin
            r_d = fit_prospect(obs, 500; version = "5b")
            means = mean(r_d)
            for sym in (:N, :Ccab, :Ccar, :Cbrown, :Cw, :Cm,
                    :ρ, :σ)
                @test isapprox(means[sym, :mean], eval(sym), rtol = 0.1)
            end
        end

        @testset "PROSPECT 5" begin
            r_d = fit_prospect(obs, 500; version = "5")
            means = mean(r_d)
            for sym in (:N, :Ccab, :Ccar, :Cw, :Cm,
                    :ρ, :σ)
                @test isapprox(means[sym, :mean], eval(sym), rtol = 0.1)
            end
        end

        @testset "PROSPECT 4" begin
            r_d = fit_prospect(obs, 500; version = "4")
            means = mean(r_d)
            for sym in (:N, :Ccab, :Cw, :Cm,
                    :ρ, :σ)
                @test isapprox(means[sym, :mean], eval(sym), rtol = 0.1)
            end
        end

    end
end

