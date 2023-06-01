# https://mathoverflow.net/questions/275831/determinant-of-correlation-matrix-of-autoregressive-model
# det_AR1(ρ, n) = (1 - ρ^2)^(n - 1)
logdet_AR1(ρ, n) = (n - 1) * log(1 - ρ^2)

# https://math.stackexchange.com/questions/975069/the-inverse-of-ar-structure-correlation-matrix-kac-murdock-szeg%C5%91-matrix
inv_AR1(ρ, n) = inv(one(ρ) - ρ^2) * SymTridiagonal([one(ρ), fill(1 + ρ^2, n-2)..., one(ρ)], fill(-ρ, n-1))

function logpdf_ar1(x, μ, σ, ρ)
    n = size(μ, 1)
    σ⁻¹ = inv(Diagonal(σ))
    s = (x - μ)
    Σ_inv = σ⁻¹ * inv_AR1(ρ, n) * σ⁻¹
    sΣs = s' * Σ_inv * s
    # Ω = Correlation; Σ = Covariance
    # det(Σ) = det(Diagonal(σ²)) * det(Ω)
    # Determinant of diagonal matrix is product of elements.
    # We do this on the variance, not the standard deviation,
    # hence σ², which is 2log(σ) in log space.
    # Product in log space --> sum.
    ldet = logdet_AR1(ρ, n) + 2 * sum(log.(σ))
    logp = n * log(2π) + ldet + sΣs
    return -0.5 * logp
end

function logpdf_ar1(x::AbstractMatrix{<:Real}, μ, σ, ρ)
    sum(map(xᵢ-> logpdf_ar1(xᵢ, μ, σ, ρ), eachcol(x)))
end

macro make_fit_prospect(name, args...)
    typedargs = [:($arg::T) for arg in args]
    namedargs = [Expr(:kw, arg, arg) for arg in args]
    priors = [:($arg ~ $(Symbol(arg, "_prior"));) for arg in args]
    return eval(quote
            function $name(obs::Spectrum, nsamples)
                opti_c = createLeafOpticalStruct(obs.λ; method = :interp)
                function myprospect($(typedargs...)) where {T}
                    leaf = LeafProspectProProperties{T}($(namedargs...))
                    _, R = prospect(leaf, opti_c)
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

N_prior = truncated(Normal(1.4, 0.2); lower = 1.0)
Ccab_prior = truncated(Normal(40, 20); lower = 0.0)
Ccar_prior = truncated(Normal(5, 20); lower = 0.0)
Canth_prior = truncated(Normal(5, 20); lower = 0.0)
Cbrown_prior = Exponential(2.0)
Cw_prior = truncated(Normal(0.01, 0.01); lower = 0.0)
Cprot_prior = truncated(Normal(0.01, 0.01); lower = 0.0)
Ccbc_prior = truncated(Normal(0.01, 0.01); lower = 0.0)
Cm_prior = truncated(Normal(0.01, 0.01); lower = 0.0)
