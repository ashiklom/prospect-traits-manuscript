using CanopyOptics

using TransformVariables, TransformedLogDensities, LogDensityProblems, LogDensityProblemsAD,
    DynamicHMC, DynamicHMC.Diagnostics, SimpleUnPack, Statistics, Random

struct ProspectProblem
    reflectance::AbstractArray{<:AbstractFloat}
end

function (problem::ProspectProblem)(θ)
    @unpack N, Cab, Cw, Cm, resid = θ
    @unpack reflectance = problem
end
