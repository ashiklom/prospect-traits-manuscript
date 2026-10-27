"""Literal transliteration of SIP_Model.m (github.com/chenminbnu/SIP-RTM-Leaf,
v1.0), including its ``eps`` guards; test oracle for ``leafrtm.rtm.leafsip``."""

import numpy as np

EPS = np.finfo(float).eps


def sip_model(data, Cab, Car, Ant, Brown, Cw, Cm):
    alpha = 600
    Kab, Kcar, Kant, KBrown, Kw, Km = data.T
    Kall = (Cab * Kab + Car * Kcar + Ant * Kant + Brown * KBrown + Cw * Kw + Cm * Km) / (Cm * alpha)
    w0 = np.exp(-Kall)
    fLMA = 2765.0 * Cm
    gLMA = 102.8 * Cm
    p = 1 - (1 - np.exp(-fLMA)) / (fLMA + EPS)
    q = 1 - 2 * np.exp(-gLMA)
    qabs = (q**2) ** 0.5
    w = w0 * (1 - p) / (1 - p * w0 + EPS)
    refl = w * (1 / 2 + q / 2 * (1 - p * w0) / (1 - qabs * p * w0))
    tran = w * (1 / 2 - q / 2 * (1 - p * w0) / (1 - qabs * p * w0))
    return refl, tran
