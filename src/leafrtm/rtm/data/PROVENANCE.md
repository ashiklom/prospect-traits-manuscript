# Provenance of the leaf RTM implementations and coefficient tables

## Code

| Module | Origin | Changes |
|---|---|---|
| `e1.py` | ashiklom/prospect-inversion-comparison @ `84504ed`, `python/prospect_np.py` and `python/prospect_pt.py` | Dropped the Swamee–Ojha approximation and the JAX kernels. The fast E1 is the NAG S13AAF Chebyshev expansion used in the original PROSPECT Fortran, via PEcAnRTM `prospect.expint.f90`. |
| `gpm.py` | the same, which ports `gpm.R`/`tav_abs` from ashiklom/rrtm (MIT) | Stokes N-layer solution split into `stokes_np`. |
| `prospect.py` | the same: a version is a coefficient table plus a constituent list (`config/models.json`) | Added PROSPECT-PRO. |
| `liberty.py` | ashiklom/rrtm `R/liberty.R` (MIT, © 2019 Alexey Shiklomanov), checked against Dawson et al. (1998) | Behavior-preserving refactors (see the module docstring), tested against a literal transliteration (`tests/reference/liberty_rrtm.py`). One correction to rrtm: the internal interface reflectance `mi` (Eq. 2) used the Fresnel reflectance for light entering the cell, not leaving it. Two variants: `liberty` (class `Liberty`) fixes a bug in the original paper, replacing Eq. 16 with the correct sum of its Eq. 15 (Eq. 16 drops a factor in the algebra, which leaves R∞ undefined for weak absorption); `liberty_dawson` (class `LibertyDawson`) keeps Eq. 16 as published. |
| `leafsip.py` | Re-implemented from Wu et al. (2021) and the reference MATLAB `SIP_Model.m` v1.0 (github.com/chenminbnu/SIP-RTM-Leaf; Y. Zeng, M. Chen, S. Wu, D. Hao; no licence stated) | Tested against a literal transliteration (`tests/reference/sip_matlab.py`). |

## Coefficient tables

| File | Contents | Source |
|---|---|---|
| `prospectd.csv` | PROSPECT-D: refractive index `n`, `talf`/`t12`/`t21`, and K for Cab, Car, Canth, Cbrown, Cw, Cm at 400–2500 nm, 1 nm | `RTM.jl/data-raw/kmat_pd.dat` and `refractive_pd.dat`, via prospect-inversion-comparison `data/coefficients/`; identical to `rrtm/data-raw/prospectd_raw.dat`. The `t21` uses the PROSPECT-D refractive index (rrtm fix `9b8a2f2`). |
| `prospect5.csv` | PROSPECT-5 and 5B: Cab, Car, Cbrown, Cw, Cm | `RTM.jl/data-raw/kmat_p5.dat` and `refractive_p45.dat` (as above) |
| `prospect4.csv` | PROSPECT-4: Cab, Cw, Cm | `RTM.jl/data-raw/kmat_p4.dat` and `refractive_p45.dat` (as above) |
| *(downloaded)* `dataSpec_PRO.txt` | PROSPECT-PRO: adds K for proteins (Cprot) and carbon-based constituents (Ccbc); the other columns are identical to PROSPECT-D | jbferet/prospect @ `9719d36`, `data-raw/dataSpec_PRO.txt` (GPL-3). It is downloaded and checksum-verified at run time (`leafrtm.rtm.prospect.pro_table`), not redistributed here. |
| `liberty.csv` | LIBERTY: K for chlorophyll, water, "albino" (dry visible), lignin+cellulose, protein at 400–2500 nm, 5 nm | `rrtm/data-raw/liberty_raw.dat` (Dawson et al. 1998 values). The model uses the first 420 rows (400–2495 nm), as rrtm does. |

`interface_terms()` checks the precomputed interface transmissivities against `tav()` on every load.

## Validation

- **PROSPECT-D** matches the independent `prosail` package (J. Gómez-Dans) to 1e-8 on five parameter sets, and matches the reference values in `rrtm/tests/testthat/test.prospect.R`.
- **PROSPECT-5/5B** match `prosail` to float32 precision; prosail stores that table in float32.
- **PROSPECT-PRO** matches `prosail` PROSPECT-D when protein and CBC absorption are passed as a combined dry-matter term.

## References

- Allen, W. A., Gausman, H. W., Richardson, A. J., Thomas, J. R. (1969). Interaction of isotropic light with a compact plant leaf. *JOSA* 59, 1376–1379.
- Allen, W. A. (1973). Transmission of isotropic light across a dielectric surface in two and three dimensions. *JOSA* 63, 664–666.
- Dawson, T. P., Curran, P. J., Plummer, S. E. (1998). LIBERTY — Modeling the effects of leaf biochemical concentration on reflectance spectra. *Remote Sens. Environ.* 65, 50–60. doi:10.1016/S0034-4257(98)00007-8
- Féret, J.-B., François, C., Asner, G. P., et al. (2008). PROSPECT-4 and 5: Advances in the leaf optical properties model separating photosynthetic pigments. *Remote Sens. Environ.* 112, 3030–3043. doi:10.1016/j.rse.2008.02.012
- Féret, J.-B., Gitelson, A. A., Noble, S. D., Jacquemoud, S. (2017). PROSPECT-D: Towards modeling leaf optical properties through a complete lifecycle. *Remote Sens. Environ.* 193, 204–215. doi:10.1016/j.rse.2017.03.004
- Féret, J.-B., Berger, K., de Boissieu, F., Malenovský, Z. (2021). PROSPECT-PRO for estimating content of nitrogen-containing leaf proteins and other carbon-based constituents. *Remote Sens. Environ.* 252, 112173. doi:10.1016/j.rse.2020.112173
- Jacquemoud, S., Verhoef, W., Baret, F., et al. (2009). PROSPECT+SAIL models: A review of use for vegetation characterization. *Remote Sens. Environ.* 113, S56–S66. doi:10.1016/j.rse.2008.01.026
- Stern, F. (1964). Transmission of isotropic radiation across an interface between two dielectrics. *Appl. Opt.* 3, 111–113.
- Stokes, G. G. (1862). On the intensity of the light reflected from or transmitted through a pile of plates. *Proc. R. Soc. Lond.* 11, 545–556.
- Wu, S., Zeng, Y., Hao, D., et al. (2021). Quantifying leaf optical properties with spectral invariants theory. *Remote Sens. Environ.* 253, 112131. doi:10.1016/j.rse.2020.112131
