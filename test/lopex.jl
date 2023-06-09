using Revise
import Pkg; Pkg.activate(".")

using ProspectTraits
using CanopyOptics
using Turing
using DataFrames, Arrow

function getobs(rn)
    metadata = DataFrame(Arrow.Table("data/ecosis-processed/lopex/metadata.arrow"))
    spectra_data = DataFrame(Arrow.Table("data/ecosis-processed/lopex/spectra.arrow"))
    observation_id = metadata[rn,:observation_id]
    observation = as_spectrum(spectra_data, observation_id)
    return observation
end

function fit_lopex(version; nsamp=1000, rn=1)
    observation = getobs(rn)
    opti_c = createLeafOpticalStruct(observation; prospect_version = version)
    turingmod = Dict(
        "4"   => ProspectTraits.prospect4_turing,
        "5"   => ProspectTraits.prospect5_turing,
        "5b"  => ProspectTraits.prospect5b_turing,
        "d"   => ProspectTraits.prospectd_turing,
        "pro" => ProspectTraits.prospectpro_turing
    )[version]
    return sample(
        turingmod(observation.values, opti_c),
        NUTS(),
        nsamp
    )
end

# Using ForwardDiff: ~190 seconds for 1000 iterations
result = fit_lopex("4")

using Zygote
Turing.setadbackend(:zygote)
rzygote = fit_lopex("4")
