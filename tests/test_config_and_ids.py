import pytest

from phylofoundry.config import validate_config
from phylofoundry.constants import DEFAULT_CONFIG
from phylofoundry.utils.id_normalization import normalize_sequence_id


def test_config_validation_good():
    cfg = {**DEFAULT_CONFIG}
    cfg["inputs"] = {**DEFAULT_CONFIG["inputs"], "faa_dir": "a", "hmm_input": "b"}
    cfg["output"] = {"outdir": "o"}
    out = validate_config(cfg)
    assert out["inputs"]["faa_dir"] == "a"


def test_config_validation_bad_type():
    cfg = {**DEFAULT_CONFIG}
    cfg["inputs"] = {**DEFAULT_CONFIG["inputs"], "faa_dir": "a", "hmm_input": "b"}
    cfg["output"] = {"outdir": "o"}
    cfg["resources"] = {"cpu": "abc"}
    with pytest.raises(SystemExit):
        validate_config(cfg)


def test_id_normalization_after_last_pipe():
    assert normalize_sequence_id("genome|prot A", mode="after_last_pipe") == "prot_A"
    assert normalize_sequence_id("a|b|c", mode="after_last_pipe") == "c"
