# import Pkg; Pkg.activate("ProspectTraits.jl")
using Revise

using ProspectTraits
using CSV
using DataFrames
using Unitful

lopex_01 = CSV.read(
    "../../spectra_db/lopex/spectra/lopex_01.csvy",
    DataFrame;
    comment = "#"
)

obs = Spectrum(
    lopex_01[:,"wavelengths"]u"nm",
    Array(select(lopex_01, Not(:wavelengths)))
    # lopex_01[:,"lopex_01.001"]
)

result = ProspectTraits.fit_prospect(obs, 5000)
