if ~isinteractive()
    using LoggingExtras
    logpath = mkpath("logs")
    logfile = "$logpath/fit-spectradb.log"
    global_logger(MinLevelLogger(FileLogger(logfile), Logging.Info))
end

using ProspectTraits

using Arrow
using DataFrames
using Base.Filesystem

using FLoops

# Read all metadata
function get_observation_ids(dataset_id)
    metadata_path = "data/ecosis-processed/$dataset_id/metadata.arrow"
    @assert isfile(metadata_path)
    df = DataFrame(Arrow.Table(metadata_path))
    return DataFrame(dataset_id = dataset_id, observation_id = df[:, :observation_id])
end
datasets_all = readdir("data/ecosis-processed")
skip_datasets = ["accp", "foster_beetle"]
datasets = setdiff(datasets_all, skip_datasets)
metadata_all = vcat([get_observation_ids(did) for did in datasets]...)

const versions = ("pro", "d", "5b", "5", "4")
for row in eachrow(metadata_all), version in versions
    dataset_id = row[:dataset_id]
    obs_id = row[:observation_id]
    try
        fit_observation(dataset_id, obs_id, version)
    catch err
        println("Error in PROSPECT inversion...")
    end
end

# dataset_id = "accp"
# observation_id = "accp|92BHIS10BW1|1992"
# spectra_df = spectra_data
