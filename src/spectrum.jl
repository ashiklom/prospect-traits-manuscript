struct Spectrum{Tλ <: AbstractArray{<:Union{Real, Quantity}},
        Tvalues <: AbstractArray{<:Real}}
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
    ""
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
