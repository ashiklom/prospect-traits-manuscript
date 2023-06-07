using Test

using ProspectTraits
using Unitful

@testset "ProspectTraits" begin
    @testset "Spectrum" begin
        @test_warn "No wavelength units given. Assuming nm." Spectrum(400.0:1.0:2500.0, rand(2101))
        @test_nowarn Spectrum((400.0:1.0:2500.0)u"nm", rand(2101))
    end
end
