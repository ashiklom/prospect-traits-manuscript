using LoggingExtras
logpath = mkpath("logs")
logfile = "$logpath/fit-lopex.log"
global_logger(MinLevelLogger(FileLogger(logfile), Logging.Info))

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
        nsamp = 500, overwrite = false, progress = false)
    outdir = mkpath("$data_basedir/results/raw/$dataset_id/prospect-$prospect_version")
    outfile = "$outdir/$observation_id.result"
    if isfile(outfile) && ~overwrite
        @info "$outfile already exists! Skipping..."
        return outfile
    end
    observation = as_spectrum(spectra_data, observation_id)
    samples = fit_prospect(observation, nsamp; version = prospect_version, progress = progress)
    serialize(outfile, samples)
    return outfile
end

const versions = ("pro", "d", "5b", "5", "4")
@info "Begin loop"
@floop for row in eachrow(metadata), version in versions
    fit_row_save(row[:observation_id], spectra_data, version; nsamp = 2000)
end
