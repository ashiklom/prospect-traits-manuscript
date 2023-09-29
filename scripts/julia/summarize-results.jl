using Revise
import Pkg; Pkg.activate(".")

using ProspectTraits
using Glob
using Serialization

using ProgressBars

resultfiles = glob("*/*/*.result", "data/results/raw")

iter = ProgressBar(resultfiles)
summaries = [summarize_fit(file) for file in iter]

params = vcat([s[1] for s in summaries]...)
corrs = vcat([s[2] for s in summaries]...)

serialize("data/results/params", params)
serialize("data/results/corrs", corrs)
