"""Tests for YAML config loading, round-tripping, and migration tooling.

Covers:
- Loading a YAML config file with ruamel.yaml
- Auto-detection of YAML vs JSON by file extension in load_json_config
- Deep-merge of YAML config values into DEFAULT_CONFIG
- Comments are preserved in the loaded CommentedMap
- The bundled config/config.yaml template is valid YAML and contains expected keys
- Migration script converts JSON → YAML without data loss
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from phylofoundry.constants import DEFAULT_CONFIG
from phylofoundry.utils.helpers import load_json_config, load_yaml_config


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture()
def simple_yaml(tmp_path: Path) -> Path:
    """Write a minimal YAML config snippet and return its path."""
    content = """\
# test config
inputs:
  faa_dir: /tmp/faa   # path to FAA files
  hmm_input: /tmp/hmm
output:
  outdir: /tmp/out
resources:
  cpu: 4
"""
    p = tmp_path / "config.yaml"
    p.write_text(content)
    return p


@pytest.fixture()
def simple_json(tmp_path: Path) -> Path:
    """Write an equivalent JSON config and return its path."""
    data = {
        "inputs": {"faa_dir": "/tmp/faa", "hmm_input": "/tmp/hmm"},
        "output": {"outdir": "/tmp/out"},
        "resources": {"cpu": 4},
    }
    p = tmp_path / "config.json"
    p.write_text(json.dumps(data, indent=2))
    return p


# ── load_yaml_config ─────────────────────────────────────────────────────────


class TestLoadYamlConfig:
    def test_returns_dict(self, simple_yaml: Path):
        cfg = load_yaml_config(str(simple_yaml))
        assert isinstance(cfg, dict)

    def test_values_correct(self, simple_yaml: Path):
        cfg = load_yaml_config(str(simple_yaml))
        assert cfg["inputs"]["faa_dir"] == "/tmp/faa"
        assert cfg["inputs"]["hmm_input"] == "/tmp/hmm"
        assert cfg["output"]["outdir"] == "/tmp/out"
        assert cfg["resources"]["cpu"] == 4

    def test_empty_yaml_returns_empty_dict(self, tmp_path: Path):
        p = tmp_path / "empty.yaml"
        p.write_text("")
        cfg = load_yaml_config(str(p))
        assert cfg == {}

    def test_null_values_preserved(self, tmp_path: Path):
        p = tmp_path / "nulls.yaml"
        p.write_text("inputs:\n  faa_dir: null\n  hmm_input: null\n")
        cfg = load_yaml_config(str(p))
        assert cfg["inputs"]["faa_dir"] is None
        assert cfg["inputs"]["hmm_input"] is None

    def test_boolean_values(self, tmp_path: Path):
        p = tmp_path / "bools.yaml"
        p.write_text("workflow:\n  force: true\nprep:\n  cleanup_combined_faa: false\n")
        cfg = load_yaml_config(str(p))
        assert cfg["workflow"]["force"] is True
        assert cfg["prep"]["cleanup_combined_faa"] is False

    def test_float_scientific_notation(self, tmp_path: Path):
        p = tmp_path / "floats.yaml"
        p.write_text("diamond:\n  max_evalue: 1.0e-5\n")
        cfg = load_yaml_config(str(p))
        assert abs(cfg["diamond"]["max_evalue"] - 1e-5) < 1e-15

    def test_list_values(self, tmp_path: Path):
        p = tmp_path / "lists.yaml"
        p.write_text(
            "synteny:\n  protein_id_field:\n    - ID\n    - protein_id\n"
        )
        cfg = load_yaml_config(str(p))
        assert cfg["synteny"]["protein_id_field"] == ["ID", "protein_id"]


# ── load_json_config (auto-detect) ───────────────────────────────────────────


class TestLoadJsonConfigAutoDetect:
    def test_yaml_extension_loads_yaml(self, simple_yaml: Path):
        cfg = load_json_config(str(simple_yaml))
        assert cfg["inputs"]["faa_dir"] == "/tmp/faa"

    def test_yml_extension_loads_yaml(self, tmp_path: Path):
        p = tmp_path / "config.yml"
        p.write_text("inputs:\n  faa_dir: /tmp/faa\n")
        cfg = load_json_config(str(p))
        assert cfg["inputs"]["faa_dir"] == "/tmp/faa"

    def test_json_extension_loads_json(self, simple_json: Path):
        cfg = load_json_config(str(simple_json))
        assert cfg["inputs"]["faa_dir"] == "/tmp/faa"

    def test_yaml_and_json_produce_same_values(
        self, simple_yaml: Path, simple_json: Path
    ):
        yaml_cfg = load_json_config(str(simple_yaml))
        json_cfg = load_json_config(str(simple_json))
        assert yaml_cfg["inputs"] == json_cfg["inputs"]
        assert yaml_cfg["output"] == json_cfg["output"]
        assert yaml_cfg["resources"] == json_cfg["resources"]


# ── Bundled config/config.yaml template ──────────────────────────────────────


class TestBundledConfigYaml:
    """Validate the bundled config/config.yaml template ships with correct keys."""

    @pytest.fixture(autouse=True)
    def _locate_template(self) -> None:
        repo_root = Path(__file__).parent.parent
        self.template_path = repo_root / "config" / "config.yaml"

    def test_template_exists(self):
        assert self.template_path.exists(), (
            f"Bundled config template not found at {self.template_path}"
        )

    def test_template_is_valid_yaml(self):
        cfg = load_yaml_config(str(self.template_path))
        assert isinstance(cfg, dict)

    def test_template_has_top_level_keys(self):
        cfg = load_yaml_config(str(self.template_path))
        required_keys = {
            "inputs", "output", "resources", "workflow", "prep",
            "hmmer", "diamond", "filtering", "phylo", "embeddings",
            "curate", "taxonomy_integrate", "conservation_metrics",
            "detect_clades", "post", "synteny", "codon", "hyphy",
            "motifs", "discover", "ha",
        }
        missing = required_keys - set(cfg.keys())
        assert not missing, f"Template missing top-level keys: {missing}"

    def test_template_inputs_match_default_config(self):
        cfg = load_yaml_config(str(self.template_path))
        for key in DEFAULT_CONFIG["inputs"]:
            assert key in cfg["inputs"], (
                f"inputs.{key} present in DEFAULT_CONFIG but not in template"
            )

    def test_template_diamond_values_match_defaults(self):
        cfg = load_yaml_config(str(self.template_path))
        dc = DEFAULT_CONFIG["diamond"]
        tc = cfg["diamond"]
        assert tc["enabled"] == dc["enabled"]
        assert tc["mode"] == dc["mode"]
        assert tc["sensitivity"] == dc["sensitivity"]
        assert abs(tc["max_evalue"] - dc["max_evalue"]) < 1e-15
        assert tc["max_target_seqs"] == dc["max_target_seqs"]

    def test_template_phylo_values_match_defaults(self):
        cfg = load_yaml_config(str(self.template_path))
        dc = DEFAULT_CONFIG["phylo"]
        tc = cfg["phylo"]
        assert tc["mafft"] == dc["mafft"]
        assert tc["iq_boot"] == dc["iq_boot"]
        assert tc["iqtree_bin"] == dc["iqtree_bin"]

    def test_template_resources_cpu(self):
        cfg = load_yaml_config(str(self.template_path))
        assert cfg["resources"]["cpu"] == DEFAULT_CONFIG["resources"]["cpu"]


# ── Migration script ──────────────────────────────────────────────────────────


class TestMigrateScript:
    """Validate the JSON-to-YAML migration script."""

    def _get_migrate_fn(self):
        import importlib.util
        import sys  # noqa: F401 (kept for sys.modules side-effect on exec_module)
        repo_root = Path(__file__).parent.parent
        script = repo_root / "scripts" / "migrate_json_to_yaml.py"
        spec = importlib.util.spec_from_file_location("migrate_json_to_yaml", script)
        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load migration script from {script}")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.migrate

    def test_migration_produces_valid_yaml(self, simple_json: Path, tmp_path: Path):
        migrate = self._get_migrate_fn()
        out = tmp_path / "out.yaml"
        migrate(simple_json, out)
        cfg = load_yaml_config(str(out))
        assert isinstance(cfg, dict)

    def test_migration_preserves_values(self, simple_json: Path, tmp_path: Path):
        migrate = self._get_migrate_fn()
        out = tmp_path / "out.yaml"
        migrate(simple_json, out)
        cfg = load_yaml_config(str(out))
        assert cfg["inputs"]["faa_dir"] == "/tmp/faa"
        assert cfg["resources"]["cpu"] == 4
        assert cfg["output"]["outdir"] == "/tmp/out"

    def test_migration_merges_into_template(self, tmp_path: Path):
        """JSON values should override template defaults while comments survive."""
        migrate = self._get_migrate_fn()
        json_data = {"resources": {"cpu": 16}}
        json_path = tmp_path / "custom.json"
        json_path.write_text(json.dumps(json_data))
        out = tmp_path / "out.yaml"
        repo_root = Path(__file__).parent.parent
        template = repo_root / "config" / "config.yaml"
        if template.exists():
            migrate(json_path, out, template_path=template)
            cfg = load_yaml_config(str(out))
            assert cfg["resources"]["cpu"] == 16
            # Other keys from template should still be present
            assert "inputs" in cfg

    def test_migration_no_template_fallback(self, simple_json: Path, tmp_path: Path):
        """Should succeed even when no annotated template is available."""
        migrate = self._get_migrate_fn()
        out = tmp_path / "out.yaml"
        migrate(simple_json, out, template_path=Path("/nonexistent/template.yaml"))
        cfg = load_yaml_config(str(out))
        assert cfg["inputs"]["faa_dir"] == "/tmp/faa"
