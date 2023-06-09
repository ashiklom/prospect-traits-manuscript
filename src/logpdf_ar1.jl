# https://mathoverflow.net/questions/275831/determinant-of-correlation-matrix-of-autoregressive-model
# det_AR1(ρ, n) = (1 - ρ^2)^(n - 1)
logdet_AR1(ρ, n) = (n - 1) * log(1 - ρ^2)

# https://math.stackexchange.com/questions/975069/the-inverse-of-ar-structure-correlation-matrix-kac-murdock-szeg%C5%91-matrix
inv_AR1(ρ, n) = inv(one(ρ) - ρ^2) * SymTridiagonal([one(ρ), fill(1 + ρ^2, n-2)..., one(ρ)], fill(-ρ, n-1))

# Simplified tridiagonal multiplication. Added here for later experiments, but
# not used right now because it's not very reliable.
function sᵀΣs(s, σ::Number, ρ)
    n = size(s, 1)
    ρ² = ρ^2
    τ = σ^(-2)
    c = 1.0 / (1.0 - ρ²)
    Ω⁻¹d = [c, fill(1 + ρ², n-2)..., c]
    Σ⁻¹d = Ω⁻¹d .* τ
    Σ⁻¹s = -ρ * τ
    return sum(Σ⁻¹d .* s.^2) + 2 * sum(Σ⁻¹s .* s[1:(n-1)] .* s[2:n])
end

function logpdf_ar1(x::AbstractVector, μ, σ::AbstractVector, ρ)
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

function logpdf_ar1(x::AbstractVector, μ, σ::Number, ρ)
    n = size(μ, 1)
    s = (x - μ)
    Σ_inv = σ^(-2) * inv_AR1(ρ, n)
    sΣs = s' * Σ_inv * s
    # Normally, sum along the diagonal. However, if σ is scalar,
    # sum downweights it...so instead, multiply by `n`.
    ldet = logdet_AR1(ρ, n) + 2n * log.(σ)
    logp = n * log(2π) + ldet + sΣs
    return -0.5 * logp
end

function logpdf_ar1(x::AbstractMatrix, μ, σ, ρ)
    sum(map(xᵢ-> logpdf_ar1(xᵢ, μ, σ, ρ), eachcol(x)))
end

