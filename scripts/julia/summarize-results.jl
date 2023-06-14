using Revise
import Pkg; Pkg.activate(".")

using ProspectTraits
using Glob

resultfiles = glob("*/*/*.result", "data/results/raw")

summaries = [summarize_fit(file) for file in resultfiles]
params = vcat([s[1] for s in summaries]...)
corrs = vcat([s[2] for s in summaries]...)
