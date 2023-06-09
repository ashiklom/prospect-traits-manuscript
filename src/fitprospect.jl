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

function prospect_defaults(version)
    default_values = LeafProspectProProperties{Float64}(Cm = 0.01)
    default_names = fieldnames(LeafProspectProProperties)
    defaults_all = (; [(name, getfield(default_values, name)) for name in default_names]...)
    prospect_params = Dict(
        "4"   => (:N, :Ccab, :Cw, :Cm),
        "5"   => (:N, :Ccab, :Ccar, :Cw, :Cm),
        "5b"  => (:N, :Ccab, :Ccar, :Cbrown, :Cw, :Cm),
        "d"   => (:N, :Ccab, :Ccar, :Canth, :Cbrown, :Cw, :Cm),
        "pro" => (:N, :Ccab, :Ccar, :Canth, :Cbrown, :Cw, :Ccbc, :Cprot)
    )
    defaults_prosp = defaults_all[prospect_params[version]]
    defaults_other = (rsd = 0.01, ρ = 0.9)
    return (; defaults_prosp..., defaults_other...)
end

const prospect_models = Dict(
    "4"   => prospect4_turing,
    "5"   => prospect5_turing,
    "5b"  => prospect5b_turing,
    "d"   => prospectd_turing,
    "pro" => prospectpro_turing
)

function fit_prospect(obs::Spectrum, nsamples::Int;
        version = "pro", sampler = NUTS(), kwargs...)
    opti_c = createLeafOpticalStruct(obs; prospect_version = version)
    turingmod = prospect_models[version]
    init_params = prospect_defaults(version)
    return sample(
        turingmod(obs.values, opti_c),
        sampler, 
        nsamples,
        init_params = init_params;
        kwargs...
    )
end

function optim_prospect(obs::Spectrum, version::String)
    opti_c = createLeafOpticalStruct(obs; prospect_version = version)
    turingmod = prospect_models[version]
    return optimize(turingmod(obs.values, opti_c), MAP())
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
