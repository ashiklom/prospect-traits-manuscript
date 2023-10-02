using Revise
import Pkg; Pkg.activate(".")

using ProspectTraits
using Glob
using CSV
using ProgressBars

resultfiles = glob("*/*/*.result", "data/results/raw")

iter = ProgressBar(resultfiles)
summaries = [summarize_fit(file) for file in iter]

params = vcat([s[1] for s in summaries]...)
corrs = vcat([s[2] for s in summaries]...)

CSV.write("data/results/params.csv", params)
CSV.write("data/results/corrs.csv", corrs)
