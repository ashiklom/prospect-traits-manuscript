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
                    resid ~ InverseGamma(1.0, 0.2)
                    mod = myprospect($(args...))
                    obs_refl ~ MvNormal(mod, resid * I)
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

# # Original implementation (for reference)
# function fit_prospect(obs::Spectrum, nsamples)
#     opti_c = createLeafOpticalStruct(obs.λ; method = :interp)
#     function prospectpro(N::T, Ccab::T, Ccar::T, Canth::T,
#             Cbrown::T, Cw::T, Cprot::T, Ccbc::T
#     ) where {T}
#         leaf = LeafProspectProProperties{T}(
#             N=N, Ccab=Ccab, Ccar=Ccar, Canth=Canth,
#             Cbrown=Cbrown, Cw=Cw, Cprot=Cprot, Ccbc=Ccbc
#         )
#         _, R = prospect(leaf, opti_c)
#         return R
#     end
#
#     @model function turingmod(obs_refl)
#         N ~ N_prior
#         Cab ~ Cab_prior
#         Car ~ Car_prior
#         Canth ~ Canth_prior
#         Cbrown ~ Cbrown_prior
#         Cw ~ Cw_prior
#         Cprot ~ Cprot_prior
#         Ccbc ~ Ccbc_prior
#         resid ~ InverseGamma(1, 0.2)
#         mod = prospect4(N, Cab, Cw, Cm)
#         obs_refl ~ MvNormal(mod, resid * I)
#     end
#
#     function sample_model(n)
#         sample(turingmod(obs.values), NUTS(), n)
#     end
#
#     # Do the real sampling --- 5000 iterations
#     samples = sample_model(nsamples)
#     return samples
# end
