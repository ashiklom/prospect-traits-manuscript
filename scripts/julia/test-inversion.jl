import Pkg; Pkg.activate(".")
using Revise

using ProspectTraits

using Arrow
using DataFrames

using Unitful

using Serialization
using Base.Filesystem

data_basedir = "data/"
dataset_id = "lopex"

const metadata = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/metadata.arrow"))
const spectra_data = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/spectra.arrow"))

function fit_row(metadata, spectra_data, rn; nsamp = 500)
    observation_id = metadata[rn, :observation_id]
    observation = as_spectrum(spectra_data, observation_id)
    return fit_prospectpro(observation, nsamp)
end

# 249.54 seconds
r1 = fit_row(metadata, spectra_data, 1; nsamp = 500)

using StatsPlots
plot(r1)
