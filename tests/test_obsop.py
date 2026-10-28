import numpy as np
import pytest

from leafrtm.obsop import Block, clean_grid, make_obsdata, operator

MODEL_WL = np.arange(400, 2501, 1.0)


def test_selection_on_model_grid_is_identity():
    m = np.random.default_rng(0).uniform(size=len(MODEL_WL))
    obs = np.arange(450, 2400, 7.0)
    idx, w = operator(MODEL_WL, obs)
    np.testing.assert_allclose(np.sum(w * m[idx], 1), m[(obs - 400).astype(int)], atol=1e-15)


def test_linear_interpolation_off_grid():
    m = 0.1 + 1e-4 * MODEL_WL  # linear in wavelength: interpolation is exact
    obs = np.array([400.4, 1000.25, 2499.9])
    idx, w = operator(MODEL_WL, obs)
    np.testing.assert_allclose(np.sum(w * m[idx], 1), 0.1 + 1e-4 * obs)


@pytest.mark.parametrize("fwhm", [2.0, 10.0, 25.0])
def test_flat_spectrum_preserved(fwhm):
    m = np.full(len(MODEL_WL), 0.37)
    idx, w = operator(MODEL_WL, np.arange(420, 2480, 10.0), fwhm=fwhm)
    np.testing.assert_allclose(np.sum(w * m[idx], 1), 0.37)
    np.testing.assert_allclose(w.sum(1), 1.0)


def test_coarse_model_grid():
    lib_wl = np.arange(400, 2496, 5.0)  # LIBERTY grid
    obs = np.arange(400, 2496, 1.0)
    idx, w = operator(lib_wl, obs)
    m = lib_wl / 1000
    np.testing.assert_allclose(np.sum(w * m[idx], 1), obs / 1000)


def test_clean_grid_sorts_and_averages_duplicates():
    wl = np.array([500, 400, 450, 450, np.nan, 600])
    v = np.array([[5, 50], [4, 40], [1, 10], [3, 30], [9, 9], [np.nan, 1]], float)
    w2, v2, seg = clean_grid(wl, v)
    np.testing.assert_array_equal(w2, [400, 450, 500])
    np.testing.assert_array_equal(v2, [[4, 40], [2, 20], [5, 50]])


def test_make_obsdata_drops_out_of_range_and_marks_starts():
    wl = np.r_[np.arange(350, 600, 1.0), np.arange(1000, 2600, 2.0)]
    b = Block("R", wl, np.linspace(0.1, 0.5, len(wl)))
    od = make_obsdata([b, Block("PA", wl, np.linspace(0.3, 1.0, len(wl)))], MODEL_WL)
    assert od.wavelength.min() >= 400 and od.wavelength.max() <= 2500
    starts = od.wavelength[od.start]
    assert list(starts) == [400.0, 400.0]  # gaps decorrelate through distance, not restarts
    assert od.quantities == ("R", "PA")
    R = np.full(len(MODEL_WL), 0.2)
    T = np.full(len(MODEL_WL), 0.3)
    pred = od.predict_np(R, T)
    np.testing.assert_allclose(pred[od.qcode == 0], 0.2)
    np.testing.assert_allclose(pred[od.qcode == 3], -np.log10(0.2))


def test_absorbance_transform():
    od = make_obsdata([Block("A", np.arange(400, 500, 1.0), np.zeros(100))], MODEL_WL)
    pred = od.predict_np(np.full(len(MODEL_WL), 0.2), np.full(len(MODEL_WL), 0.3))
    np.testing.assert_allclose(pred, 0.5)


def test_usable_bands_trims_noisy_detector_ends_only():
    """Noise-dominated ends (and impossible values) are dropped; the interior and
    smooth ends are kept, even where the interior has sharp features."""
    from leafrtm.prep import usable_bands

    rng = np.random.default_rng(0)
    wl = np.arange(400.0, 2501.0)
    y = 0.3 + 0.2 * np.tanh((wl - 720) / 20) + rng.normal(0, 1e-4, wl.size)  # smooth, with a red edge
    assert usable_bands(y, wl, "R").all()
    noisy = y.copy()
    noisy[wl > 2440] += rng.normal(0, 0.15, (wl > 2440).sum())  # wild SWIR tail
    noisy[wl < 420] += rng.normal(0, 0.05, (wl < 420).sum())    # noisy blue end
    noisy[1200] = -0.4                                           # impossible value inside
    k = usable_bands(noisy, wl, "R")
    assert not k[wl > 2445].any() and not k[wl < 415].any()
    assert not k[1200] and k[(wl > 450) & (wl < 2350) & (wl != wl[1200])].all()
