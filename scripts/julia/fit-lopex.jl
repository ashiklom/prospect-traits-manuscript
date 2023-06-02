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

function fit_row_save(metadata, spectra_data, rn;
        nsamp = 500, overwrite = false)
    observation_id = metadata[rn, :observation_id]
    outdir = mkpath("$data_basedir/results/raw/$dataset_id/prospect-pro")
    outfile = "$outdir/$observation_id.result"
    if isfile(outfile) && ~overwrite
        println("$outfile already exists! Skipping...")
        return outfile
    end
    observation = as_spectrum(spectra_data, observation_id)
    samples = fit_prospectpro(observation, nsamp)
    serialize(outfile, samples)
    return outfile
end

# Try the first row...
fit_row_save(metadata, spectra_data, 1; nsamp = 2000)
