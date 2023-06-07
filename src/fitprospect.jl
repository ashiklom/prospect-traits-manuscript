# https://mathoverflow.net/questions/275831/determinant-of-correlation-matrix-of-autoregressive-model
# det_AR1(ρ, n) = (1 - ρ^2)^(n - 1)
logdet_AR1(ρ, n) = (n - 1) * log(1 - ρ^2)

# https://math.stackexchange.com/questions/975069/the-inverse-of-ar-structure-correlation-matrix-kac-murdock-szeg%C5%91-matrix
inv_AR1(ρ, n) = inv(one(ρ) - ρ^2) * SymTridiagonal([one(ρ), fill(1 + ρ^2, n-2)..., one(ρ)], fill(-ρ, n-1))

function logpdf_ar1(x::AbstractVector, μ, σ, ρ)
    n = size(μ, 1)
    s = (x - μ)
    σ⁻¹ = inv(Diagonal(σ))
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

function logpdf_ar1(x::AbstractMatrix, μ, σ, ρ)
    sum(map(xᵢ-> logpdf_ar1(xᵢ, μ, σ, ρ), eachcol(x)))
end

@model function prospectpro_turing(obs_refl, opti_c)
    N ~ N_prior
    Ccab ~ Ccab_prior
    Ccar ~ Ccar_prior
    Canth ~ Canth_prior
    Cbrown ~ Cbrown_prior
    Cw ~ Cw_prior
    Ccbc ~ Ccbc_prior
    Cprot ~ Cprot_prior
    # σ_a ~ σa_prior
    # σ_b ~ σb_prior
    σ ~ σ_prior
    ρ ~ ρ_prior
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Canth=Canth, Cbrown=Cbrown,
        Cw=Cw, Ccbc=Ccbc, Cprot=Cprot,
        Cm=zero(N)
    )
    _, pred = prospect(leaf, opti_c)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    σᵥ = fill(σ, size(pred, 1))
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σᵥ, ρ)
end

@model function prospectd_turing(obs_refl, opti_c)
    N ~ N_prior
    Ccab ~ Ccab_prior
    Ccar ~ Ccar_prior
    Canth ~ Canth_prior
    Cbrown ~ Cbrown_prior
    Cw ~ Cw_prior
    Cm ~ Cm_prior
    # σ_a ~ σa_prior
    # σ_b ~ σb_prior
    σ ~ σ_prior
    ρ ~ ρ_prior
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Canth=Canth, Cbrown=Cbrown,
        Cw=Cw, Cm=Cm,
        Ccbc=zero(N), Cprot=zero(N)
    )
    _, pred = prospect(leaf, opti_c)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    σᵥ = fill(σ, size(pred, 1))
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σᵥ, ρ)
end

@model function prospect5b_turing(obs_refl, opti_c)
    N ~ N_prior
    Ccab ~ Ccab_prior
    Ccar ~ Ccar_prior
    Cbrown ~ Cbrown_prior
    Cw ~ Cw_prior
    Cm ~ Cm_prior
    # σ_a ~ σa_prior
    # σ_b ~ σb_prior
    σ ~ σ_prior
    ρ ~ ρ_prior
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Cbrown=Cbrown,
        Cw=Cw, Cm=Cm,
        Canth = zero(N), Ccbc=zero(N), Cprot=zero(N)
    )
    _, pred = prospect(leaf, opti_c)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    σᵥ = fill(σ, size(pred, 1))
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σᵥ, ρ)
end

@model function prospect5_turing(obs_refl, opti_c)
    N ~ N_prior
    Ccab ~ Ccab_prior
    Ccar ~ Ccar_prior
    Cw ~ Cw_prior
    Cm ~ Cm_prior
    # σ_a ~ σa_prior
    # σ_b ~ σb_prior
    σ ~ σ_prior
    ρ ~ ρ_prior
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Cw=Cw, Cm=Cm,
        Cbrown=zero(N), Canth=zero(N), Ccbc=zero(N), Cprot=zero(N)
    )
    _, pred = prospect(leaf, opti_c)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    σᵥ = fill(σ, size(pred, 1))
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σᵥ, ρ)
end

@model function prospect4_turing(obs_refl, opti_c)
    N ~ N_prior
    Ccab ~ Ccab_prior
    Cw ~ Cw_prior
    Cm ~ Cm_prior
    # σ_a ~ σa_prior
    # σ_b ~ σb_prior
    σ ~ σ_prior
    ρ ~ ρ_prior
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Cw=Cw, Cm=Cm,
        Ccar=zero(N), Cbrown=zero(N), Canth=zero(N), Ccbc=zero(N), Cprot=zero(N)
    )
    _, pred = prospect(leaf, opti_c)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    σᵥ = fill(σ, size(pred, 1))
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σᵥ, ρ)
end

function fit_prospect(obs::Spectrum, nsamples::Int;
        version = "pro", sampler = NUTS(), kwargs...)
    dλ = diff(obs.λ)
    λ_windows = vcat(obs.λ[1] - dλ[1], obs.λ[1:(end-1)] .+ dλ, obs.λ[end] + dλ[end])
    if version == "5b"
        opti_c = createLeafOpticalStruct(λ_windows; prospect_version = "5")
    else
        opti_c = createLeafOpticalStruct(λ_windows; prospect_version = version)
    end
    turingmod = Dict(
        "4"   => prospect4_turing,
        "5"   => prospect5_turing,
        "5b"  => prospect5b_turing,
        "d"   => prospectd_turing,
        "pro" => prospectpro_turing
    )[version]
    return sample(
        turingmod(obs.values, opti_c),
        sampler, 
        nsamples; 
        kwargs...
    )
end

N_prior = truncated(Normal(1.4, 0.2); lower = 1.0)
Ccab_prior = truncated(Normal(40, 20); lower = 0.0)
Ccar_prior = truncated(Normal(5, 20); lower = 0.0)
Canth_prior = truncated(Normal(5, 20); lower = 0.0)
Cbrown_prior = Exponential(2.0)
Cw_prior = truncated(Normal(0.01, 0.01); lower = 0.0)
Cprot_prior = truncated(Normal(0.01, 0.01); lower = 0.0)
Ccbc_prior = truncated(Normal(0.01, 0.01); lower = 0.0)
Cm_prior = truncated(Normal(0.01, 0.01); lower = 0.0)

σ_prior = Exponential(0.1)
# σa_prior = Exponential(0.005)
# σb_prior = Exponential(0.02)
ρ_prior = Beta(12.0, 1.1)
