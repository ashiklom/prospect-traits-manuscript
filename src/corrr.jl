"""
Convert correlation matrix (or correlation matrix computed from MCMC chains)
to a tidy dataframe, a la R's `corrr::stretch`.
"""
function stretch(x::ChainDataFrame)
    xdf = DataFrame(x)
    pars = String.(xdf[:,:parameters])
    mat = Matrix(xdf[:,2:end])
    return stretch(mat, pars)
end

function stretch(x::Matrix, pars)
    par1 = vcat([repeat([pars[i]], size(pars,1)-i) for i in eachindex(pars)]...)
    par2 = vcat([pars[(i+1):end] for i in eachindex(pars)]...)
    vecL = vcat([x[(i+1):end,i] for i in 1:size(x,1)]...)
    return DataFrame(var1 = par1, var2 = par2, corr = vecL)
end

function stretch(x::Matrix)
    pars = ["var$i" for i in 1:size(x,1)]
    return stretch(x, pars)
end

# x = cor(result)
# xm = Matrix(DataFrame(x)[:,2:end])
# xdf = stretch(cor(result))
# xm2 = unstretch(xdf)
# xm2 == xm
# xm == stretch()

"""
Convert tidy correlation dataframe back to a correlation matrix.
"""
function unstretch(xdf::DataFrame)
    @assert unique(xdf[:,:var1])[2:end] == unique(xdf[:,:var2])[1:(end-1)]
    pars = unique(vcat(xdf[:,:var1]..., xdf[:,:var2]...))
    M = zeros(length(pars), length(pars))
    M[diagind(M)] .= 1.0
    for i in eachindex(pars)
        p = pars[i]
        y = subset(xdf, :var1 => x -> x .== p)[:,:corr]
        M[(i+1):end,i] = y
        M[i,(i+1):end] = y
    end
    return M
end
