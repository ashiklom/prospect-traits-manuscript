function fit_prospect(obs::Spectrum, nsamples)
    opti_c = createLeafOpticalStruct(obs.λ; method = :interp)
    function prospect4(N::T, Cab::T, Cw::T, Cm::T) where {T}
        leaf = LeafProspectProProperties{T}(N=N, Ccab=Cab, Cw=Cw, Cm=Cm)
        _, R = prospect(leaf, opti_c)
        return R
    end

    @model function turingmod(obs_refl)
        N ~ truncated(Normal(1.4, 0.2); lower = 1.0, upper = 4.0)
        Cab ~ truncated(Normal(40, 20); lower = 0.0, upper = 120.0)
        Cw ~ truncated(Normal(0.01, 0.01); lower = 0.0, upper = 0.1)
        Cm ~ truncated(Normal(0.01, 0.01); lower = 0.0, upper = 0.1)
        resid ~ InverseGamma(1, 0.2)
        mod = prospect4(N, Cab, Cw, Cm)
        obs_refl ~ MvNormal(mod, resid * I)
    end

    function sample_model(n)
        sample(turingmod(obs.values), NUTS(), n)
    end

    # Do the real sampling --- 5000 iterations
    samples = sample_model(nsamples)
    return samples
end
