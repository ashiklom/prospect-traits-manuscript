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
    pred = prospectpro(opti_c, N, Ccab, Ccar, Canth, Cbrown, Cw,
        Ccbc, Cprot)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σ, ρ)
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
    pred = prospectd(opti_c, N, Ccab, Ccar, Canth, Cbrown, Cw, Cm)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σ, ρ)
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
    pred = prospect5b(opti_c, N, Ccab, Ccar, Cbrown, Cw, Cm)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σ, ρ)
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
    pred = prospect5(opti_c, N, Ccab, Ccar, Cw, Cm)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σ, ρ)
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
    pred = prospect4(opti_c, N, Ccab, Cw, Cm)
    # Heteroskedastic variance model 
    # σ = σ_a .* pred .+ σ_b
    Turing.@addlogprob! logpdf_ar1(obs_refl, pred, σ, ρ)
end

function fit_prospect(obs::Spectrum, nsamples::Int;
        version = "pro", sampler = NUTS(), kwargs...)
    opti_c = createLeafOpticalStruct(obs; prospect_version = version)
    prospect_models = Dict(
        "4"   => prospect4_turing,
        "5"   => prospect5_turing,
        "5b"  => prospect5b_turing,
        "d"   => prospectd_turing,
        "pro" => prospectpro_turing
    )
    turingmod = prospect_models[version]
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
