import Pkg; Pkg.activate(".")
using Revise

using ProspectTraits

using Arrow
using DataFrames

using Unitful

using Serialization
using Base.Filesystem

const data_basedir = "../data/"
const dataset_id = "lopex"

const metadata = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/metadata.arrow"))
const spectra_data = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/spectra.arrow"))

function fit_row(metadata, spectra_data, rn;
        nsamp = 500, overwrite = false)
    observation_id = metadata[rn, :observation_id]
    outdir = mkpath("$data_basedir/results/raw/$dataset_id/")
    outfile = "$outdir/$observation_id.result"
    if isfile(outfile) && ~overwrite
        println("$outfile already exists! Skipping...")
        return outfile
    end
    observation = as_spectrum(spectra_data, observation_id)
    samples = fit_prospect(observation, nsamp)
    serialize(outfile, samples)
    return outfile
end

r1 = fit_row(metadata, spectra_data, 1; nsamp = 500)
r1b = fit_row(metadata, spectra_data, 1; nsamp = 500)
r2 = fit_row(metadata, spectra_data, 2; nsamp = 500)
