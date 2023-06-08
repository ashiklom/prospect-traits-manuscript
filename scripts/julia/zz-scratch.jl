import Pkg; Pkg.activate(".")
using Revise

using ProspectTraits
using Unitful
using CanopyOptics

const waves = (400.0:1.0:2500.0) * u"nm"
const opti_c = createLeafOpticalStruct(waves; method = :interp)
leaf = LeafProspectProProperties{Float64}
prosp(leaf) = prospect(leaf, opti_c)
refl = prosp(leaf(Ccab = 40.0))[2] + (randn(2101) * 0.003)

const obs_c = Spectrum(waves, refl)
m = ProspectTraits.prospectpro_turing(obs_c.values, opti_c)

using Turing
samps = sample(m, NUTS(), 500)

# init_params = optimize(m, MAP())
# samps = sample(m, NUTS(100, 0.65), 500, init_params = init_params.values.array)
samps = sample(m, NUTS(), 500; save_state = true)  # 131 sec

using StatsPlots
plot(samps)

refl_3 = hcat([prosp(leaf(Ccab = 35.0))[2] + (randn(2101) * 0.003) for _ in 1:3]...)
const obs_c3 = Spectrum(waves, refl_3)
m3 = ProspectTraits.prospectpro_turing(obs_c3.values, opti_c)
samps = sample(m, NUTS())
# init_params = optimize(m3, MAP())

using DynamicHMC
samps3 = sample(m, DynamicNUTS(), 500; save_state = true)

fit_prospectpro(obs_c3, 500) # 377.39

################################################################################
observation

################################################################################

using CanopyOptics
using Unitful

using GLMakie

const waves = (400.0:1.0:2500.0) * u"nm"
const opti_c = createLeafOpticalStruct(waves; method = :interp)

function plotmat(x, y)
    fig = Figure()
    ax = Axis(fig[1,1])
    λ = ustrip(waves)
    for i in eachindex(x)
        lines!(ax, λ, y[:,i])
    end
    return fig
end

leaf = LeafProspectProProperties{Float64}
prosp(leaf) = prospect(leaf, opti_c)
results = map(x -> prosp(leaf(N=x))[2], range(1.0, 10.0, 50))
# N can be anything

results = map(x -> prosp(leaf(Ccab=x))[2], range(0.0, 400, 50))
# Cab can be anything

xx = range(0.0, 400, 50)
results = hcat(map(x -> prosp(leaf(Ccar=x))[2], xx)...)
# plotmat(xx, results)
# Car can be anything

xx = range(0.0, 500, 500)
results = hcat(map(x -> prosp(leaf(Canth=x))[2], xx)...)
# Canth can be anything

xx = range(0.0, 500, 500)
results = hcat(map(x -> prosp(leaf(Cbrown=x))[2], xx)...)
# Cbrown can be anything

xx = range(0.0, 5.0, 500)
results = hcat(map(x -> prosp(leaf(Cw=x))[2], xx)...)
# Cw can be anything

xx = range(0.0, 5.0, 500)
results = hcat(map(x -> prosp(leaf(Cm=x))[2], xx)...)

xx = range(0.0, 5.0, 500)
results = hcat(map(x -> prosp(leaf(Ccbc=x))[2], xx)...)

xx = range(0.0, 5.0, 500)
leaf(N=1.4)
println("Leaf: ", leaf(N=1.4))
results = hcat(map(x -> prosp(leaf(Cprot=x))[2], xx)...)

################################################################################
leaf(; params...)

################################################################################

using Distributions
using LinearAlgebra

obs = [
    1.0 1.1 1.2 1.3
    2.0 2.1 2.2 2.3
    3.0 3.1 3.2 3.3
]

dist = MvNormal(zeros(3), diagm(ones(3)))

logpdf(dist, obs)

################################################################################

using Plots
using Distributions

Nprior = truncated(Normal(1.4, 0.2); lower = 1.0)
plot(0.0:0.01:2.5, x -> pdf(Nprior, x))

plot(0.0:0.1:10, x -> pdf(Cauchy(0, 0.1), x))

################################################################################
Base.@kwdef struct mystruct
    a = 3
    b = 4
    c = 5
end

macro mymacro(args...)
    namedargs = [Expr(:kw, arg, arg) for arg in args]
    println("Arguments: $namedargs")
    return eval(quote
            function makestruct($(args...))
                mystruct($(namedargs...))
            end
        end)
end
@mymacro a
makestruct(Nothing)
# makestruct(Nothing)

c = Nothing
eval(:(mystruct($(:c)=$c)))

ex = :(mystruct(b=2))
dump(:(mystruct(b=2)))

exprs = (Expr(:kw, :b, Nothing), Expr(:kw, :c, Nothing))
eval(Expr(:call, :mystruct, exprs...))
eval(:(mystruct($(exprs...))))


mystruct
Meta.show_sexpr(mystruct(b=2, c=1))

myvals = (a=1, c=2)
mystruct((b=2, c=2, a=3)...).a

function myfunc(; a=3, b=6)
    return a - b
end

################################################################################
using CanopyOptics
using Statistics
using StatsBase
using Distributions
using LinearAlgebra

const opti_c = createLeafOpticalStruct(400.0:1.0:2500.0; method = :interp)
const leaf_c = LeafProspectProProperties{Float64}()
T, R = prospect(leaf_c, opti_c)

ρ = autocor(R, [1])

MvNormal(R, (ρ * ρ') * I)

autocor(R)

################################################################################
[0:4 -1:3]
using LinearAlgebra

H(n::Int) = abs.(hcat([collect((1-i):(n-i)) for i in 1:n]...))
σ = collect(1.1:1.0:5.5)
ρ = 0.7
Ω = ρ .^ H(5)

Diagonal(σ^2) * Ω

################################################################################

# function fit_row_save(metadata, spectra_data, rn;
#         nsamp = 500, overwrite = false)
#     observation_id = metadata[rn, :observation_id]
#     outdir = mkpath("$data_basedir/results/raw/$dataset_id/")
#     outfile = "$outdir/$observation_id.result"
#     if isfile(outfile) && ~overwrite
#         println("$outfile already exists! Skipping...")
#         return outfile
#     end
#     observation = as_spectrum(spectra_data, observation_id)
#     samples = fit_prospectpro(observation, nsamp)
#     serialize(outfile, samples)
#     return outfile
# end

# using LinearAlgebra
# using PDMats
# mod = randn(10)
# σ = PDiagMat(abs.(5 .* mod .+ 3))
# ρ = 0.7
# H = ProspectTraits.Hmat(10)
# Ω = PDMat(ρ .^ H)
# Σ = X_A_Xt(σ, Ω)

# 78.28 seconds
# r1out = deserialize(r1)

# r1b = fit_row(metadata, spectra_data, 1; nsamp = 500)
# r2 = fit_row(metadata, spectra_data, 2; nsamp = 500)

################################################################################

macro make_fit_prospect(name, args...)
    typedargs = [:($arg::T) for arg in args]
    namedargs = [Expr(:kw, arg, arg) for arg in args]
    priors = [:($arg ~ $(Symbol(arg, "_prior"));) for arg in args]
    return eval(quote
            function $name(obs::Spectrum, nsamples)
                opti_c = createLeafOpticalStruct(obs.λ; method = :interp)
                function myprospect($(typedargs...)) where {T}
                    leaf = LeafProspectProProperties{T}($(namedargs...))
                    _, R = try 
                        prospect(leaf, opti_c)
                    catch
                        println("Error in PROSPECT evaluation. Current leaf:")
                        println(leaf)
                        nothing, fill(999.9, size(opti_c.r12, 1))
                    end
                    return R
                end
                @model function turingmod(obs_refl)
                    $(priors...)
                    σ_a ~ Exponential(0.005)
                    σ_b ~ Exponential(0.02)
                    ρ ~ Beta(12, 1.1)
                    pred = myprospect($(args...))
                    # Heteroskedastic variance model 
                    σ = σ_a .* pred .+ σ_b
                    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σ, ρ)
                end
                return sample(turingmod(obs.values), NUTS(), nsamples)
            end
            export $name
        end)
end

@make_fit_prospect fit_prospect4   N Ccab Cw Cm
@make_fit_prospect fit_prospect5   N Ccab Ccar Cw Cm
@make_fit_prospect fit_prospect5b  N Ccab Ccar Cbrown Cw Cm
@make_fit_prospect fit_prospectD   N Ccab Ccar Canth Cbrown Cw Cm
@make_fit_prospect fit_prospectpro N Ccab Ccar Canth Cbrown Cw Cprot Ccbc

function sᵀΣs(s, σ, ρ)
    n = size(σ, 1)
    ρ² = ρ^2
    c = 1.0 / (1.0 - ρ²)
    Ω⁻¹d = [c, fill(1 + ρ², n-2)..., c]
    Σ⁻¹d = Ω⁻¹d .* (σ.^2)
    Σ⁻¹s = -ρ .* σ[1:(n-1)] .* σ[2:n]
    return sum(Σ⁻¹d .* s.^2) + 2 * sum(Σ⁻¹s .* s[1:(n-1)] .* s[2:n])
end

################################################################################

outfile = fit_row_save(metadata[2, :observation_id], spectra_data, "pro";
    nsamp = 2000, progress = true, overwrite = true)
deserialize(outfile)

outfile = fit_row_save(metadata[3, :observation_id], spectra_data, "pro";
    nsamp = 2000, progress = true, overwrite = true)
