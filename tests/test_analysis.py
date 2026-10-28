"""Suspect-group flagging and exclusions in the trait analyses."""

import numpy as np
import pandas as pd

from leafrtm.analysis.compute import excluded, exclusion_rules, suspect_groups


def _leaf(rng):
    """Two runs, one trait; dataset c (2012) is estimated 2.7x too high, the rest within ~10 %."""
    rows = []
    for run in ("r1", "r2"):
        for ds, year, factor in (("a", 2011, 1.0), ("b", 2015, 1.1), ("c", 2011, 1.0), ("c", 2012, 2.7),
                                 ("d", 2018, 0.95)):
            for i in range(6):
                m = rng.uniform(20, 60)
                e = m * factor * np.exp(rng.normal(0, 0.05))
                rows.append({"observation_id": f"{ds}:{year}:{i}", "run_id": run, "trait": "chl_area",
                             "relative": False, "q50": e, "q05": e * 0.9, "q95": e * 1.1, "measured": m,
                             "covered90": bool(e * 0.9 <= m <= e * 1.1), "dataset_id": ds})
    leaf = pd.DataFrame(rows)
    ctx = leaf[["observation_id", "dataset_id"]].drop_duplicates().assign(
        year=lambda d: d["observation_id"].str.split(":").str[1].astype(int))
    return leaf, ctx


def test_flags_only_the_group_that_misses_consistently():
    leaf, ctx = _leaf(np.random.default_rng(1))
    s = suspect_groups(leaf, ctx)
    flagged = s[s["flagged"]]
    assert list(zip(flagged["dataset_id"], flagged["year"])) == [("c", 2012)]
    assert flagged["deviation_ratio"].iat[0] > 2


def test_exclusions_by_dataset_year_trait_and_scope():
    leaf, ctx = _leaf(np.random.default_rng(1))
    s = suspect_groups(leaf, ctx)
    cfg = {"exclude": [{"dataset": "a", "year": 2011, "traits": ["chl_area"], "reason": "units", "scope": "all"}]}
    rules = exclusion_rules(cfg, s)
    assert len(rules) == 1
    m_all = excluded(leaf, ctx, rules, scopes=("all",))
    assert set(leaf.loc[m_all, "observation_id"].str[:6]) == {"a:2011"}
    assert not excluded(leaf, ctx, rules, scopes=("ecology",)).any()
    auto = exclusion_rules({**cfg, "auto_exclude_flagged": True}, s)
    eco = excluded(leaf, ctx, auto)  # default scopes: ecology and all
    assert set(leaf.loc[eco, "observation_id"].str[:6]) == {"a:2011", "c:2012"}
