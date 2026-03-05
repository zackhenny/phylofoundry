import numpy as np
from phylofoundry.utils.ha import write_ha_outputs


def test_ha_outputs_created(tmp_path):
    seqs = {"a|p1": "ACDEFGHIK", "b|p2": "ACDEFGHIK"}
    scores = {k: np.linspace(0.1, 1.0, len(v)) for k, v in seqs.items()}
    masks = {k: (np.arange(len(v)) % 2).astype("uint8") for k, v in seqs.items()}
    attention = tmp_path / "attention"
    summary = tmp_path / "summary"
    write_ha_outputs(str(attention), str(summary), "HMM1", seqs, scores, masks, 1, 3, {"call_mode": "percentile", "percentile": 0.2})
    assert (attention / "HMM1.ha_sites.tsv").exists()
    assert (summary / "HMM1.ha_counts.tsv").exists()
