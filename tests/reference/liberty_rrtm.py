"""Literal transliteration of ``liberty()`` from ashiklom/rrtm (R/liberty.R),
used as the test oracle for ``leafrtm.rtm.liberty``. ``n_iter`` defaults to
rrtm's 50 fixed-point iterations; pass a large value to get the converged
solution. ``me_override`` and ``mi_override`` (arrays over wavelength) replace
rrtm's interface reflectances when given (``leafrtm.rtm.liberty`` corrects rrtm's ``mi``)."""

import numpy as np


def liberty_rrtm(K, D, xu, thick, baseline, element, c_factor, w_factor, l_factor,
                 p_factor, n_iter=50, me_override=None, mi_override=None):
    k_chloro, k_water, ke, k_ligcell, k_protein = K.T
    nwl = 420
    reflR = np.zeros(nwl)
    transR = np.zeros(nwl)
    RR = np.zeros(nwl)
    for i in range(nwl):
        coeff = D * (baseline + k_chloro[i] * c_factor + k_water[i] * w_factor
                     + ke[i] * element + k_ligcell[i] * l_factor + k_protein[i] * p_factor)
        N1 = 1.4891 - 0.0005 * i
        N0 = 1.0
        width = np.pi / 180
        alpha = np.arange(1, 91) * np.pi / 180
        beta = np.arcsin(N0 / N1 * np.sin(alpha))
        plus = alpha + beta
        dif = alpha - beta
        refl = 0.5 * (np.sin(dif) ** 2 / np.sin(plus) ** 2 + np.tan(dif) ** 2 / np.tan(plus) ** 2)
        me = 2 * np.sum(refl * np.sin(alpha) * np.cos(alpha) * width)
        critical = np.arcsin(N0 / N1) * 180 / np.pi
        jj = np.arange(1, critical + 1e-12)  # R: seq(1, critical)
        jj = jj[jj <= critical]
        alpha = jj * np.pi / 180
        beta = np.arcsin(N0 / N1 * np.sin(alpha))
        plus = alpha + beta
        dif = alpha - beta
        refl = 0.5 * (np.sin(dif) ** 2 / np.sin(plus) ** 2 + np.tan(dif) ** 2 / np.tan(plus) ** 2)
        mint = np.sum(refl * np.sin(alpha) * np.cos(alpha) * width)
        mi = (1 - np.sin(critical * np.pi / 180) ** 2) + 2 * mint
        if me_override is not None:
            me, mi = me_override[i], mi_override[i]
        M = 2 * (1 - (coeff + 1) * np.exp(-coeff)) / coeff**2
        Tt = ((1 - mi) * M) / (1 - (mi * M))
        x = xu / (1 - (1 - 2 * xu) * Tt)
        a = me * Tt + x * Tt - me - Tt - x * me * Tt
        b = 1 + x * me * Tt - 2 * x**2 * me**2 * Tt
        c = 2 * me * x**2 * Tt - x * Tt - 2 * x * me
        R = 0.5
        for _ in range(n_iter):
            R = -(a * R**2 + c) / b
        rb = (2 * x * me) + (x * Tt) - (x * Tt * 2 * x * me)
        tb = np.sqrt(((R - rb) * (1 - (R * rb))) / R)
        whole = np.trunc(thick)
        fraction = thick - whole
        top = tb ** (1 + fraction) * ((((1 + tb) ** 2) - rb**2) ** (1 - fraction))
        bot1 = (1 + tb) ** (2 * (1 - fraction)) - rb**2
        bot2 = 1 + (64 / 3) * fraction * (fraction - 0.5) * (fraction - 1) * 0.001
        tif = top / (bot1 * bot2)
        rif = (1 + rb**2 - tb**2 - np.sqrt((1 + rb**2 - tb**2) ** 2 - 4 * rb**2 * (1 - tif**2))) / (2 * rb)
        if whole >= 2:
            prev_t, prev_r = 1.0, 0.0
            for _ in range(int(whole) - 1):
                cur_t = (prev_t * tb) / (1 - (prev_r * rb))
                cur_r = prev_r + (((prev_t * prev_t) * rb) / (1 - (prev_r * rb)))
                prev_t, prev_r = cur_t, cur_r
        else:
            cur_t, cur_r = 1.0, 0.0
        transR[i] = (cur_t * tif) / (1 - (rif * cur_r))
        reflR[i] = cur_r + ((cur_t**2 * rif) / (1 - (rif * cur_r)))
        RR[i] = R
    return reflR, transR, RR
