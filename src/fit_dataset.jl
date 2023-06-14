function fit_observation(dataset_id, observation_id, prospect_version;
        outdir_base = "data/results/raw", nsamp = 2000,
        overwrite = false, progress = false)
    outdir = mkpath("$outdir_base/$dataset_id/prospect-$prospect_version")
    outfile = "$outdir/$observation_id.result"
    logtag = "[$observation_id, $prospect_version]"
    if isfile(outfile) && ~overwrite
        @info "$logtag Skipping b/c already complete"
        return outfile
    end
    @info "$logtag Beginning inversion"
    specfile = "data/ecosis-processed/$dataset_id/spectra.arrow"
    observation = as_spectrum(specfile, observation_id)
    samples = fit_prospect(observation, nsamp; version = prospect_version, progress = progress)
    @info "$logtag Inversion complete! Saving results"
    serialize(outfile, samples)
    return outfile
end
