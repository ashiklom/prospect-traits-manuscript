# using LoggingExtras
# logpath = mkpath("logs")
# logfile = "$logpath/fit-lopex.log"
# global_logger(MinLevelLogger(FileLogger(logfile), Logging.Info))

using ProspectTraits

using Arrow
using DataFrames
using JSON

using Unitful
using Optim

using FLoops

using Serialization
using Base.Filesystem

data_basedir = "data/"
dataset_id = "lopex"

const metadata = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/metadata.arrow"))
const spectra_data = DataFrame(Arrow.Table("$data_basedir/ecosis-processed/$dataset_id/spectra.arrow"))

# observation_id = metadata[1, :observation_id]
# prospect_version = "pro"

function optim_row_save(observation_id, spectra_data, prospect_version;
        nsamp = 500, overwrite = false, progress = false)
    outdir = mkpath("$data_basedir/results/raw-optim/$dataset_id/prospect-$prospect_version")
    outfile = "$outdir/$observation_id.json"
    if isfile(outfile) && ~overwrite
        @info "$outfile already exists! Skipping..."
        return outfile
    end
    observation = as_spectrum(spectra_data, observation_id)
    result = optim_prospect(observation, prospect_version)
    results_dict = Dict(result.values)
    output = Dict(
        "dataset_id" => dataset_id,
        "observation_id" => observation_id,
        "estimates" => results_dict,
        "diagonstics" => Dict(
            "lp" => result.lp,
            "converged" => Optim.converged(result.optim_result)
        )
    )
    open(outfile, "w") do f
        write(f, json(output))
    end
    return outfile
end

const versions = ("pro", "d", "5b", "5", "4")
@info "Begin loop"
@floop for row in eachrow(metadata), version in versions
    fit_row_save(row[:observation_id], spectra_data, version; nsamp = 2000)
end
