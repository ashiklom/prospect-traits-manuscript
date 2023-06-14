# fname = "data/results/raw/angers/prospect-4/an03r0001.result"
# fname = "data/results/raw/angers/prospect-4/an03r0002.result"
# summarize_fit(fname)

function summarize_fit(fname)
    rx = r"(?:data/)?(?:results/)?(?:raw/)?" *
      r"(?<dataset_id>.*?)/" * 
      r"prospect-(?<version>.*?)/" *
      r"(?<observation_id>.*)\.result$"
    m = match(rx, fname)
    dataset_id = String(m[:dataset_id])
    version = String(m[:version])
    result = deserialize(fname)
    param_summary = DataFrame(summarystats(result))
    param_quant = DataFrame(quantile(result; q = [0.025, 0.975]))
    params = outerjoin(param_summary, param_quant, on = :parameters)
    params[!,:dataset_id] .= dataset_id
    params[!,:version] .= version
    select!(params, :dataset_id, :version, Not([:dataset_id, :version]))
    corrdf = stretch(cor(result))
    corrdf[!,:dataset_id] .= dataset_id
    corrdf[!,:version] .= version
    select!(corrdf, :dataset_id, :version, Not([:dataset_id, :version]))
    return params, corrdf
end
