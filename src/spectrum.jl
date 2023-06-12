struct Spectrum{Tλ <: AbstractArray{typeof(400.0u"nm")}, Tvalues}
    λ::Tλ
    values::Tvalues
end

function Spectrum(λ::AbstractArray{<:Real}, values)
    @warn "No wavelength units given. Assuming nm."
    Spectrum(λ*u"nm", values)
end

"""
    as_spectrum(spectra_df, observation_id)

Convert spectra DataFrame row to Spectrum object suitable for inversion

# Arguments
- `spectra_df::DataFrame`: Data frame containing spectra
- `observation_id::String`: Observation ID

"""
function as_spectrum(spectra_df, observation_id)
    spec_obs = subset(
        spectra_df,
        :observation_id => x -> x .== observation_id,
        :spectral_measurement => x -> x .== "reflectance"
    )[:, [:wavelength_nm, :spectra_id, :value]]
    spec_wide = unstack(spec_obs, :spectra_id, :value)
    waves = Array{Float64}(spec_wide[:, :wavelength_nm])*u"nm"
    values = Array{Float64}(spec_wide[:, Not(:wavelength_nm)])
    Spectrum(waves, values)
end

function CanopyOptics.createLeafOpticalStruct(obs::Spectrum; prospect_version = "pro")
    dλ = diff(obs.λ)
    λ_windows = vcat(obs.λ[1] - dλ[1], obs.λ[1:(end-1)] .+ dλ, obs.λ[end] + dλ[end])
    v = prospect_version == "5b" ? "5" : prospect_version
    return createLeafOpticalStruct(λ_windows; prospect_version = v)
end
