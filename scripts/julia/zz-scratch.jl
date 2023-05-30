using Distributions
using LinearAlgebra

obs = [
    1.0 1.1 1.2 1.3
    2.0 2.1 2.2 2.3
    3.0 3.1 3.2 3.3
]

dist = MvNormal(zeros(3), diagm(ones(3)))

logpdf(dist, obs)

################################################################################

using Plots
using Distributions

Nprior = truncated(Normal(1.4, 0.2); lower = 1.0)
plot(0.0:0.01:2.5, x -> pdf(Nprior, x))

plot(0.0:0.1:10, x -> pdf(Cauchy(0, 0.1), x))

################################################################################
Base.@kwdef struct mystruct
    a = 3
    b = 4
    c = 5
end

macro mymacro(args...)
    namedargs = [Expr(:kw, arg, arg) for arg in args]
    println("Arguments: $namedargs")
    return eval(quote
            function makestruct($(args...))
                mystruct($(namedargs...))
            end
        end)
end
@mymacro a
makestruct(Nothing)
# makestruct(Nothing)

c = Nothing
eval(:(mystruct($(:c)=$c)))

ex = :(mystruct(b=2))
dump(:(mystruct(b=2)))

exprs = (Expr(:kw, :b, Nothing), Expr(:kw, :c, Nothing))
eval(Expr(:call, :mystruct, exprs...))
eval(:(mystruct($(exprs...))))


mystruct
Meta.show_sexpr(mystruct(b=2, c=1))

myvals = (a=1, c=2)
mystruct((b=2, c=2, a=3)...).a

function myfunc(; a=3, b=6)
    return a - b
end

################################################################################
using CanopyOptics
using Statistics
using StatsBase
using Distributions
using LinearAlgebra

const opti_c = createLeafOpticalStruct(400.0:1.0:2500.0; method = :interp)
const leaf_c = LeafProspectProProperties{Float64}()
T, R = prospect(leaf_c, opti_c)

ρ = autocor(R, [1])

MvNormal(R, (ρ * ρ') * I)

autocor(R)
