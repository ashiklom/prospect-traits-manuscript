module ProspectTraits

using CanopyOptics
using Turing
using Optim
using Distributions
using Unitful
using LinearAlgebra
using DataFrames

include("spectrum.jl")
include("logpdf_ar1.jl")
include("prospect.jl")
include("fitprospect.jl")

export Spectrum,
       as_spectrum,
       createLeafOpticalStruct,
       fit_prospect,
       optim_prospect

end # module ProspectTraits
