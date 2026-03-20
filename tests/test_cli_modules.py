"""Tests for the new modular CLI subcommand interface.

Verifies that:
- Each pipeline step can be invoked as a standalone subcommand.
- Step-specific arguments are applied to the resolved config.
- Optional steps are auto-enabled when invoked as a subcommand.
- Utility subcommands (list-steps, plan, validate, doctor, dump-config) work.
- The legacy flag-style invocation continues to work (backward compatibility).
"""

from __future__ import annotations

import argparse
import json
import sys
from copy import deepcopy
from unittest.mock import patch

import pytest

from phylofoundry.constants import DEFAULT_CONFIG, STEPS
from phylofoundry.main import (
    _ALL_SUBCMDS,
    _STEP_ENABLE_MAP,
    _STEP_SUBCMD_MAP,
    _apply_step_args,
    _build_deps,
    _build_parser,
)


# ── Helpers ───────────────────────────────────────────────────────────────────


def _run_main(argv: list[str]) -> int:
    """Run main() with patched sys.argv; return the SystemExit code."""
    from phylofoundry.main import main

    with patch("sys.argv", ["phylofoundry"] + argv):
        with pytest.raises(SystemExit) as exc_info:
            main()
    return exc_info.value.code


def _valid_cfg(**overrides) -> dict:
    cfg = deepcopy(DEFAULT_CONFIG)
    cfg["inputs"]["faa_dir"] = "/tmp/faa"
    cfg["inputs"]["hmm_input"] = "/tmp/hmm"
    cfg["output"]["outdir"] = "/tmp/out"
    cfg.update(overrides)
    return cfg


def _parse(argv: list[str]) -> argparse.Namespace:
    """Parse *argv* through the full argument parser."""
    from phylofoundry.main import _ALL_SUBCMDS, _build_parser

    # Apply the same legacy-routing logic as main()
    if not argv or argv[0] not in _ALL_SUBCMDS:
        argv = ["run"] + argv
    ap = _build_parser()
    return ap.parse_args(argv)


# ── _ALL_SUBCMDS completeness ─────────────────────────────────────────────────


class TestSubcmdRegistry:
    def test_all_steps_have_subcommand(self):
        """Every internal step name must map to a CLI subcommand."""
        step_to_subcmd = {v: k for k, v in _STEP_SUBCMD_MAP.items()}
        for step in STEPS:
            assert step in step_to_subcmd, (
                f"Step '{step}' has no CLI subcommand mapping in _STEP_SUBCMD_MAP"
            )

    def test_all_subcmds_registered_in_all_subcmds(self):
        for cmd in _STEP_SUBCMD_MAP:
            assert cmd in _ALL_SUBCMDS

    def test_utility_subcmds_in_all_subcmds(self):
        for cmd in ("run", "list-steps", "plan", "validate", "doctor", "dump-config"):
            assert cmd in _ALL_SUBCMDS


# ── Backward compatibility (legacy flag routing) ───────────────────────────────


class TestLegacyFlagBackwardCompat:
    """Old-style ``phylofoundry --flag`` invocations must still work."""

    def test_list_steps_legacy_flag(self, capsys):
        code = _run_main(["--list-steps"])
        assert code == 0
        assert "prep" in capsys.readouterr().out

    def test_plan_legacy_flag(self, capsys):
        code = _run_main([
            "--plan",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code == 0
        assert "Execution plan" in capsys.readouterr().out

    def test_validate_config_legacy_flag_valid(self, capsys):
        code = _run_main([
            "--validate-config",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code == 0

    def test_validate_config_legacy_flag_invalid(self, capsys):
        code = _run_main(["--validate-config"])
        assert code != 0

    def test_doctor_legacy_flag(self, capsys):
        code = _run_main([
            "--doctor",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code in (0, 1)

    def test_dump_default_config_legacy_flag(self, capsys):
        code = _run_main(["--dump_default_config"])
        assert code == 0
        out = capsys.readouterr().out
        parsed = json.loads(out)
        assert "inputs" in parsed


# ── New utility subcommands ────────────────────────────────────────────────────


class TestUtilitySubcommands:
    def test_list_steps_subcommand(self, capsys):
        code = _run_main(["list-steps"])
        assert code == 0
        out = capsys.readouterr().out
        for step in STEPS:
            assert step in out

    def test_dump_config_subcommand(self, capsys):
        code = _run_main(["dump-config"])
        assert code == 0
        parsed = json.loads(capsys.readouterr().out)
        assert "inputs" in parsed
        assert "output" in parsed

    def test_plan_subcommand(self, capsys):
        code = _run_main([
            "plan",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code == 0
        assert "Execution plan" in capsys.readouterr().out

    def test_validate_subcommand_valid(self, capsys):
        code = _run_main([
            "validate",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code == 0

    def test_validate_subcommand_invalid(self, capsys):
        code = _run_main(["validate"])
        assert code != 0

    def test_doctor_subcommand(self):
        code = _run_main([
            "doctor",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code in (0, 1)


# ── Step subcommands: routing and start_at/stop_after ────────────────────────


class TestStepSubcommandRouting:
    """Verify that step subcommands set start_at and stop_after correctly."""

    @pytest.mark.parametrize("subcmd,expected_step", list(_STEP_SUBCMD_MAP.items()))
    def test_step_subcommand_sets_start_stop(self, subcmd, expected_step):
        args = _parse([subcmd])
        from phylofoundry.main import _STEP_SUBCMD_MAP as MAP

        step_internal = MAP.get(args.subcommand)
        assert step_internal == expected_step

    def test_embed_subcommand_detected(self):
        args = _parse(["embed"])
        assert args.subcommand == "embed"

    def test_phylo_subcommand_detected(self):
        args = _parse(["phylo"])
        assert args.subcommand == "phylo"

    def test_hmmer_subcommand_detected(self):
        args = _parse(["hmmer"])
        assert args.subcommand == "hmmer"

    def test_detect_clades_hyphen_subcommand(self):
        args = _parse(["detect-clades"])
        assert args.subcommand == "detect-clades"
        from phylofoundry.main import _STEP_SUBCMD_MAP
        assert _STEP_SUBCMD_MAP["detect-clades"] == "detect_clades"

    def test_score_motifs_hyphen_subcommand(self):
        args = _parse(["score-motifs"])
        assert args.subcommand == "score-motifs"

    def test_discover_motifs_hyphen_subcommand(self):
        args = _parse(["discover-motifs"])
        assert args.subcommand == "discover-motifs"

    def test_taxonomy_maps_to_taxonomy_integrate(self):
        args = _parse(["taxonomy"])
        from phylofoundry.main import _STEP_SUBCMD_MAP
        assert _STEP_SUBCMD_MAP["taxonomy"] == "taxonomy_integrate"

    def test_conservation_maps_to_conservation_metrics(self):
        args = _parse(["conservation"])
        from phylofoundry.main import _STEP_SUBCMD_MAP
        assert _STEP_SUBCMD_MAP["conservation"] == "conservation_metrics"


# ── Step-specific argument parsing ────────────────────────────────────────────


class TestStepSpecificArgs:
    """Verify that step-specific CLI args are applied to the config."""

    def test_embed_model_arg(self):
        cfg = _valid_cfg()
        args = _parse(["embed", "--outdir", "/tmp/out", "--model", "esm2_t6_8M_UR50D"])
        _apply_step_args(args, cfg)
        assert cfg["embeddings"]["model"] == "esm2_t6_8M_UR50D"

    def test_embed_device_arg(self):
        cfg = _valid_cfg()
        args = _parse(["embed", "--outdir", "/tmp/out", "--device", "cpu"])
        _apply_step_args(args, cfg)
        assert cfg["embeddings"]["device"] == "cpu"

    def test_embed_batch_size_arg(self):
        cfg = _valid_cfg()
        args = _parse(["embed", "--outdir", "/tmp/out", "--batch_size", "16"])
        _apply_step_args(args, cfg)
        assert cfg["embeddings"]["batch_size"] == 16

    def test_embed_backend_arg(self):
        cfg = _valid_cfg()
        args = _parse(["embed", "--outdir", "/tmp/out", "--backend", "transformers"])
        _apply_step_args(args, cfg)
        assert cfg["embeddings"]["backend"] == "transformers"

    def test_phylo_mafft_arg(self):
        cfg = _valid_cfg()
        args = _parse(["phylo", "--outdir", "/tmp/out", "--mafft"])
        _apply_step_args(args, cfg)
        assert cfg["phylo"]["mafft"] is True

    def test_phylo_combined_arg(self):
        cfg = _valid_cfg()
        args = _parse(["phylo", "--outdir", "/tmp/out", "--combined"])
        _apply_step_args(args, cfg)
        assert cfg["phylo"]["combined_tree"] is True

    def test_phylo_iqtree_bin_arg(self):
        cfg = _valid_cfg()
        args = _parse(["phylo", "--outdir", "/tmp/out", "--iqtree_bin", "iqtree2"])
        _apply_step_args(args, cfg)
        assert cfg["phylo"]["iqtree_bin"] == "iqtree2"

    def test_phylo_iq_boot_arg(self):
        cfg = _valid_cfg()
        args = _parse(["phylo", "--outdir", "/tmp/out", "--iq_boot", "500"])
        _apply_step_args(args, cfg)
        assert cfg["phylo"]["iq_boot"] == 500

    def test_hmmer_diamond_mode_arg(self):
        cfg = _valid_cfg()
        args = _parse(["hmmer", "--diamond_mode"])
        _apply_step_args(args, cfg)
        assert cfg["diamond"]["enabled"] is True

    def test_hmmer_diamond_query_arg(self):
        cfg = _valid_cfg()
        args = _parse(["hmmer", "--diamond_query", "/tmp/query.faa"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["diamond_query"] == "/tmp/query.faa"

    def test_taxonomy_gtdb_dir_arg(self):
        cfg = _valid_cfg()
        args = _parse(["taxonomy", "--gtdb_dir", "/tmp/gtdb"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["gtdb_dir"] == "/tmp/gtdb"

    def test_taxonomy_taxonomy_file_arg(self):
        cfg = _valid_cfg()
        args = _parse(["taxonomy", "--taxonomy_file", "/tmp/tax.tsv"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["taxonomy_file"] == "/tmp/tax.tsv"

    def test_conservation_clades_tsv_arg(self):
        cfg = _valid_cfg()
        args = _parse(["conservation", "--clades_tsv", "/tmp/clades.tsv"])
        _apply_step_args(args, cfg)
        assert cfg["conservation_metrics"]["clades_tsv"] == "/tmp/clades.tsv"

    def test_detect_clades_detect_method_arg(self):
        cfg = _valid_cfg()
        args = _parse(["detect-clades", "--detect_method", "treecluster"])
        _apply_step_args(args, cfg)
        assert cfg["detect_clades"]["detect_method"] == "treecluster"

    def test_synteny_gbk_dir_arg(self):
        cfg = _valid_cfg()
        args = _parse(["synteny", "--gbk_dir", "/tmp/gbk"])
        _apply_step_args(args, cfg)
        assert cfg["synteny"]["gbk_dir"] == "/tmp/gbk"

    def test_synteny_gff_dir_arg(self):
        cfg = _valid_cfg()
        args = _parse(["synteny", "--gff_dir", "/tmp/gff"])
        _apply_step_args(args, cfg)
        assert cfg["synteny"]["gff_dir"] == "/tmp/gff"

    def test_codon_pal2nal_cmd_arg(self):
        cfg = _valid_cfg()
        args = _parse(["codon", "--pal2nal_cmd", "/opt/pal2nal.pl"])
        _apply_step_args(args, cfg)
        assert cfg["codon"]["pal2nal_cmd"] == "/opt/pal2nal.pl"

    def test_hyphy_hyphy_bin_arg(self):
        cfg = _valid_cfg()
        args = _parse(["hyphy", "--hyphy_bin", "/usr/local/bin/hyphy"])
        _apply_step_args(args, cfg)
        assert cfg["hyphy"]["hyphy_bin"] == "/usr/local/bin/hyphy"

    def test_hyphy_hyphy_tests_arg(self):
        cfg = _valid_cfg()
        args = _parse(["hyphy", "--hyphy_tests", "RELAX,aBSREL"])
        _apply_step_args(args, cfg)
        assert cfg["hyphy"]["hyphy_tests"] == "RELAX,aBSREL"

    def test_score_motifs_motifs_arg(self):
        cfg = _valid_cfg()
        args = _parse(["score-motifs", "--motifs", "HPEVY,HPEVF"])
        _apply_step_args(args, cfg)
        assert cfg["motifs"]["enabled"] is True
        assert cfg["motifs"]["motif_list"] == ["HPEVY", "HPEVF"]

    def test_run_combined_arg(self):
        cfg = _valid_cfg()
        args = _parse(["run", "--combined"])
        _apply_step_args(args, cfg)
        assert cfg["phylo"]["combined_tree"] is True

    def test_run_motifs_arg(self):
        cfg = _valid_cfg()
        args = _parse(["run", "--motifs", "HPEVY"])
        _apply_step_args(args, cfg)
        assert cfg["motifs"]["enabled"] is True
        assert "HPEVY" in cfg["motifs"]["motif_list"]


# ── Auto-enable optional steps ────────────────────────────────────────────────


class TestAutoEnable:
    """Optional steps must be auto-enabled when invoked as a subcommand."""

    @pytest.mark.parametrize("step_internal,spec", list(_STEP_ENABLE_MAP.items()))
    def test_optional_step_is_auto_enabled(self, step_internal, spec):
        section, key = spec
        cfg = _valid_cfg()
        assert cfg.get(section, {}).get(key, False) is False, (
            f"Step '{step_internal}' should default to disabled"
        )
        cfg.setdefault(section, {})[key] = True
        assert cfg[section][key] is True

    def test_embed_auto_enabled_end_to_end(self):
        """Smoke-test that embed subcommand sets embeddings.enabled=True in config."""
        from phylofoundry.main import _STEP_ENABLE_MAP
        section, key = _STEP_ENABLE_MAP["embed"]
        cfg = _valid_cfg()
        cfg[section][key] = True  # simulate auto-enable
        assert cfg["embeddings"]["enabled"] is True

    def test_mandatory_steps_not_in_enable_map(self):
        """prep, hmmer, extract, phylo are mandatory and must NOT be in _STEP_ENABLE_MAP."""
        for step in ("prep", "hmmer", "extract", "phylo"):
            assert step not in _STEP_ENABLE_MAP, (
                f"Mandatory step '{step}' should not require auto-enabling"
            )


# ── _build_deps ────────────────────────────────────────────────────────────────


class TestBuildDeps:
    def test_prep_has_no_deps(self):
        cfg = _valid_cfg()
        assert _build_deps(cfg, "prep") == []

    def test_extract_has_no_deps(self):
        cfg = _valid_cfg()
        assert _build_deps(cfg, "extract") == []

    def test_embed_has_no_deps(self):
        cfg = _valid_cfg()
        assert _build_deps(cfg, "embed") == []

    def test_hmmer_returns_hmmer_tools(self):
        cfg = _valid_cfg()
        deps = _build_deps(cfg, "hmmer")
        assert "hmmscan" in deps
        assert "hmmsearch" in deps

    def test_hmmer_diamond_mode_returns_diamond(self):
        cfg = _valid_cfg()
        cfg["diamond"]["enabled"] = True
        deps = _build_deps(cfg, "hmmer")
        assert "diamond" in deps
        assert "hmmscan" not in deps

    def test_phylo_returns_iqtree_mafft_clipkit(self):
        cfg = _valid_cfg()
        deps = _build_deps(cfg, "phylo")
        assert "mafft" in deps
        assert "clipkit" in deps
        assert any("iqtree" in d for d in deps)

    def test_codon_returns_pal2nal(self):
        cfg = _valid_cfg()
        deps = _build_deps(cfg, "codon")
        assert any("pal2nal" in d for d in deps)

    def test_hyphy_returns_hyphy(self):
        cfg = _valid_cfg()
        deps = _build_deps(cfg, "hyphy")
        assert any("hyphy" in d for d in deps)

    def test_full_pipeline_includes_hmmer_and_phylo(self):
        cfg = _valid_cfg()
        deps = _build_deps(cfg, None)
        assert "hmmscan" in deps or "diamond" in deps
        assert any("iqtree" in d for d in deps)

    def test_curate_has_no_deps(self):
        cfg = _valid_cfg()
        assert _build_deps(cfg, "curate") == []

    def test_taxonomy_integrate_has_no_deps(self):
        cfg = _valid_cfg()
        assert _build_deps(cfg, "taxonomy_integrate") == []

    def test_conservation_metrics_has_no_deps(self):
        cfg = _valid_cfg()
        assert _build_deps(cfg, "conservation_metrics") == []

    def test_detect_clades_has_no_deps(self):
        cfg = _valid_cfg()
        assert _build_deps(cfg, "detect_clades") == []


# ── Common args on every step subcommand ─────────────────────────────────────


class TestCommonArgsOnAllSubcmds:
    """Each step subcommand must accept --config, --outdir, --cpu, --force."""

    @pytest.mark.parametrize("subcmd", list(_STEP_SUBCMD_MAP.keys()))
    def test_common_args_accepted(self, subcmd):
        args = _parse([
            subcmd,
            "--outdir", "/tmp/out",
            "--cpu", "4",
            "--force",
        ])
        assert args.outdir == "/tmp/out"
        assert args.cpu == 4
        assert args.force is True

    @pytest.mark.parametrize("subcmd", list(_STEP_SUBCMD_MAP.keys()))
    def test_config_arg_accepted(self, subcmd):
        args = _parse([subcmd, "--config", "/tmp/config.json"])
        assert args.config == "/tmp/config.json"


# ── plan subcommand respects start_at/stop_after ─────────────────────────────


class TestPlanSubcommandOptions:
    def test_plan_start_at_stop_after(self, capsys):
        code = _run_main([
            "plan",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
            "--start_at", "hmmer",
            "--stop_after", "extract",
        ])
        assert code == 0
        out = capsys.readouterr().out
        assert "Execution plan" in out


# ── run subcommand: full-pipeline args ────────────────────────────────────────


class TestRunSubcommand:
    def test_run_start_at_stop_after(self):
        args = _parse([
            "run",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
            "--start_at", "embed",
            "--stop_after", "phylo",
        ])
        assert args.start_at == "embed"
        assert args.stop_after == "phylo"

    def test_run_list_steps_flag(self, capsys):
        code = _run_main(["run", "--list-steps"])
        assert code == 0
        assert "prep" in capsys.readouterr().out

    def test_run_dump_default_config(self, capsys):
        code = _run_main(["run", "--dump_default_config"])
        assert code == 0
        parsed = json.loads(capsys.readouterr().out)
        assert "inputs" in parsed
