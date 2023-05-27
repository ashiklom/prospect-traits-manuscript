module ProspectTraits
export Spectrum, fit_prospect

using CanopyOptics
using Turing
using Distributions
using Unitful
using LinearAlgebra
using DataFrames

include("spectrum.jl")
include("fitprospect.jl")

end # module ProspectTraits
