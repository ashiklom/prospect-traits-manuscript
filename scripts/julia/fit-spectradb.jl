using Distributed
addprocs()

using Dates
using Arrow
using DataFrames
using Base.Filesystem

# Read all metadata
function get_observation_ids(dataset_id)
    metadata_path = "data/ecosis-processed/$dataset_id/metadata.arrow"
    @assert isfile(metadata_path)
    df = DataFrame(Arrow.Table(metadata_path))
    return DataFrame(dataset_id = dataset_id, observation_id = df[:, :observation_id])
end
datasets_all = readdir("data/ecosis-processed")
skip_datasets = ["accp", "foster_beetle", "wisc-leaf-trait-vine"]
datasets = setdiff(datasets_all, skip_datasets)
metadata_all = vcat([get_observation_ids(did) for did in datasets]...)

@everywhere using ProspectTraits

logpath = mkpath("logs")
dstamp = Dates.format(now(), "yyyymmdd-HHMMSS")
logfile = "$logpath/fit-spectradb-$dstamp.log"
@everywhere using LoggingExtras
@everywhere global_logger(MinLevelLogger(FileLogger($logfile; append=true), Logging.Info))

@everywhere function fit_safe(dataset_id, obs_id, version)
    result = try 
        fit_observation(dataset_id, obs_id, version)
    catch err
        if isa(err, InterruptException)
            rethrow(err)
        else
            @info "[$dataset_id, $obs_id] $err\n\nSkipping..."
            nothing
        end
    end
    return result
end

const versions = ("pro", "d", "5b", "5", "4")
items = [(row[:dataset_id], row[:observation_id], version) for row in eachrow(metadata_all) for version in versions]
@sync @distributed for item = items
    fit_safe(item...)
end
