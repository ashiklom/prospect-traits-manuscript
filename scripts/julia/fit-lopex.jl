# using Revise
@info "Activatintg package"
import Pkg; Pkg.activate(".")

@info "Loading ProspectTraits"
using ProspectTraits

@info "Loading other libraries"
using Arrow
using DataFrames

using Unitful

using FLoops

using Serialization
using Base.Filesystem

@info "Begin code"
data_basedir = "data/"
dataset_id = "lopex"

const metadata = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/metadata.arrow"))
const spectra_data = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/spectra.arrow"))

function fit_row_save(observation_id, spectra_data, prospect_version;
        nsamp = 500, overwrite = false)
    outdir = mkpath("$data_basedir/results/raw/$dataset_id/prospect-$prospect_version")
    outfile = "$outdir/$observation_id.result"
    if isfile(outfile) && ~overwrite
        println("$outfile already exists! Skipping...")
        return outfile
    end
    observation = as_spectrum(spectra_data, observation_id)
    samples = fit_prospect(observation, nsamp; version = prospect_version, progress = false)
    serialize(outfile, samples)
    return outfile
end

# Try the first row...
# @time fit_row_save(metadata[1,:observation_id], spectra_data; nsamp = 2000)

const versions = ("pro", "d", "5b", "5", "4")
# Threads.@threads for row in eachrow(metadata), version in versions
# Threads.@threads for (row, version) in Iterators.product(eachrow(metadata), versions)
@info "Begin loop"
@floop for row in eachrow(metadata), version in versions
    fit_row_save(row[:observation_id], spectra_data, version; nsamp = 2000)
end
