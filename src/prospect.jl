function prospect4(opti_c, N::T, Ccab::T, Cw::T, Cm::T) where {T}
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Cw=Cw, Cm=Cm,
        Ccar=zero(T), Cbrown=zero(T), Canth=zero(T), Ccbc=zero(T), Cprot=zero(T)
    )
    _, pred = prospect(leaf, opti_c)
    return pred
end

function prospect5(opti_c, N::T, Ccab::T, Ccar::T, Cw::T, Cm::T) where {T}
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Cw=Cw, Cm=Cm,
        Cbrown=zero(T), Canth=zero(T), Ccbc=zero(T), Cprot=zero(T)
    )
    _, pred = prospect(leaf, opti_c)
    return pred
end

function prospect5b(opti_c, N::T, Ccab::T, Ccar::T, Cbrown::T, Cw::T, Cm::T) where {T}
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Cbrown=Cbrown, Cw=Cw, Cm=Cm,
         Canth=zero(T), Ccbc=zero(T), Cprot=zero(T)
    )
    _, pred = prospect(leaf, opti_c)
    return pred
end

function prospectd(opti_c, N::T, Ccab::T, Ccar::T, Canth::T, Cbrown::T, Cw::T, Cm::T) where {T}
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Canth=Canth, Cbrown=Cbrown, Cw=Cw, Cm=Cm,
         Ccbc=zero(T), Cprot=zero(T)
    )
    _, pred = prospect(leaf, opti_c)
    return pred
end

function prospectpro(opti_c, N::T, Ccab::T, Ccar::T, Canth::T,
        Cbrown::T, Cw::T, Cprot::T, Ccbc::T) where {T}
    leaf = LeafProspectProProperties(
        N=N, Ccab=Ccab, Ccar=Ccar, Canth=Canth, Cbrown=Cbrown, Cw=Cw, 
        Cprot=Cprot, Ccbc=Ccbc,
        Cm=zero(T)
    )
    _, pred = prospect(leaf, opti_c)
    return pred
end
