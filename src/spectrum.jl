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
- `spectra_df`: `DataFrame` containing spectra, or `String` path to spectra
  dataset in Arrow format.
- `observation_id`: Observation ID

"""
function as_spectrum(spectra_df, observation_id)
    spec_obs = subset(
        spectra_df,
        :observation_id => x -> x .== observation_id,
        :spectral_measurement => x -> x .== "reflectance",
        :wavelength_nm => x -> x .>= 400.0,
        :wavelength_nm => x -> x .<= 2500.0,
        skipmissing = true
    )[:, [:wavelength_nm, :spectra_id, :value]]
    dropmissing!(spec_obs)
    @assert nrow(spec_obs) > 0 "No reflectance observations found for $observation_id"
    spec_wide = unstack(spec_obs, :spectra_id, :value)
    waves = (spec_wide[:, :wavelength_nm])*u"nm"
    values = spec_wide[:, Not(:wavelength_nm)]
    Spectrum(waves, values)
end

function as_spectrum(spectra_df::String, observation_id)
    @assert isfile(spectra_df) "$spectra_df not found"
    return as_spectrum(DataFrame(Arrow.Table(spectra_df)), observation_id)
end

function CanopyOptics.createLeafOpticalStruct(obs::Spectrum; prospect_version = "pro")
    dλ = diff(obs.λ)
    λ_windows = vcat(obs.λ[1] - dλ[1], obs.λ[1:(end-1)] .+ dλ, obs.λ[end] + dλ[end])
    v = prospect_version == "5b" ? "5" : prospect_version
    return createLeafOpticalStruct(λ_windows; prospect_version = v)
end
