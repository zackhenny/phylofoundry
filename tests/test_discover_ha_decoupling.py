import pandas as pd
import pytest

from phylofoundry.tasks.discover import discover_motifs


def test_discover_requires_ha_when_enabled(tmp_path):
    cfg = {"discover": {"enabled": True, "use_ha": True}}
    with pytest.raises(SystemExit, match="HA artifacts missing"):
        discover_motifs(cfg, str(tmp_path), str(tmp_path), hmm_keep=None, force=True)


def test_discover_runs_without_ha_when_disabled(tmp_path):
    cfg = {"discover": {"enabled": True, "use_ha": False}}
    out = discover_motifs(cfg, str(tmp_path), str(tmp_path), hmm_keep=None, force=True)
    assert isinstance(out, pd.DataFrame)
    assert out.empty
