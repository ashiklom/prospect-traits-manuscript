using Distributions
using LinearAlgebra

obs = [
    1.0 1.1 1.2 1.3
    2.0 2.1 2.2 2.3
    3.0 3.1 3.2 3.3
]

dist = MvNormal(zeros(3), diagm(ones(3)))

logpdf(dist, obs)

