"""Tests for the preflight / planning / doctor / validate-config modes."""

import argparse
import sys
from copy import deepcopy
from unittest.mock import patch

import pytest

from phylofoundry.constants import DEFAULT_CONFIG, STEPS
from phylofoundry.preflight import (
    print_list_steps,
    print_plan,
    run_doctor,
    validate_config_only,
)


# ── Helpers ───────────────────────────────────────────────────────────────────


def _valid_cfg(**overrides) -> dict:
    """Return a full valid config dict (deep copy of DEFAULT_CONFIG + required fields)."""
    cfg = deepcopy(DEFAULT_CONFIG)
    cfg["inputs"]["faa_dir"] = "/tmp/faa"
    cfg["inputs"]["hmm_input"] = "/tmp/hmm"
    cfg["output"]["outdir"] = "/tmp/out"
    cfg.update(overrides)
    return cfg


# ── print_list_steps ──────────────────────────────────────────────────────────


class TestListSteps:
    def test_outputs_all_steps(self, capsys):
        print_list_steps()
        out = capsys.readouterr().out
        for step in STEPS:
            assert step in out

    def test_includes_description(self, capsys):
        print_list_steps()
        out = capsys.readouterr().out
        # The prep step's description should appear
        assert "Combine input FAA files" in out

    def test_shows_optional_tag_for_optional_steps(self, capsys):
        print_list_steps()
        out = capsys.readouterr().out
        # embed is optional
        assert "[optional]" in out

    def test_shows_depends_on(self, capsys):
        print_list_steps()
        out = capsys.readouterr().out
        # hmmer depends on prep
        assert "prep" in out


# ── print_plan ────────────────────────────────────────────────────────────────


class TestPrintPlan:
    def test_smoke(self, capsys):
        cfg = _valid_cfg()
        print_plan(cfg)
        out = capsys.readouterr().out
        assert "Execution plan" in out

    def test_enabled_steps_shown(self, capsys):
        cfg = _valid_cfg()
        print_plan(cfg)
        out = capsys.readouterr().out
        # prep, hmmer, extract, phylo should always be enabled
        assert "prep" in out
        assert "hmmer" in out
        assert "extract" in out
        assert "phylo" in out

    def test_disabled_optional_steps_shown(self, capsys):
        cfg = _valid_cfg()
        # embeddings is optional and off by default
        assert cfg["embeddings"]["enabled"] is False
        print_plan(cfg)
        out = capsys.readouterr().out
        assert "embed" in out

    def test_start_at_stop_after_respected(self, capsys):
        cfg = _valid_cfg()
        cfg["workflow"]["start_at"] = "hmmer"
        cfg["workflow"]["stop_after"] = "extract"
        print_plan(cfg)
        out = capsys.readouterr().out
        # prep comes before hmmer so it should be DISABLED
        assert "DISABLED" in out
        # hmmer and extract are within the range so ENABLED should appear
        assert "ENABLED" in out

    def test_tip_shown_when_disabled_optional_steps_exist(self, capsys):
        cfg = _valid_cfg()
        print_plan(cfg)
        out = capsys.readouterr().out
        # The tip message references enabling via config flag
        assert "enabled" in out and "config" in out.lower()


# ── validate_config_only ──────────────────────────────────────────────────────


class TestValidateConfigOnly:
    def test_valid_config_returns_true(self, capsys):
        cfg = _valid_cfg()
        result = validate_config_only(cfg)
        assert result is True
        out = capsys.readouterr().out
        assert "valid" in out.lower()

    def test_missing_faa_dir_returns_false(self, capsys):
        cfg = _valid_cfg()
        cfg["inputs"]["faa_dir"] = None
        result = validate_config_only(cfg)
        assert result is False
        out = capsys.readouterr().out
        assert "faa_dir" in out

    def test_missing_hmm_input_returns_false(self, capsys):
        cfg = _valid_cfg()
        cfg["inputs"]["hmm_input"] = None
        result = validate_config_only(cfg)
        assert result is False
        out = capsys.readouterr().out
        assert "hmm_input" in out

    def test_missing_outdir_returns_false(self, capsys):
        cfg = _valid_cfg()
        cfg["output"]["outdir"] = None
        result = validate_config_only(cfg)
        assert result is False
        out = capsys.readouterr().out
        assert "outdir" in out

    def test_invalid_start_at_returns_false(self, capsys):
        cfg = _valid_cfg()
        cfg["workflow"]["start_at"] = "nonexistent_step"
        result = validate_config_only(cfg)
        assert result is False
        out = capsys.readouterr().out
        assert "start_at" in out

    def test_invalid_stop_after_returns_false(self, capsys):
        cfg = _valid_cfg()
        cfg["workflow"]["stop_after"] = "nonexistent_step"
        result = validate_config_only(cfg)
        assert result is False

    def test_start_after_stop_returns_false(self, capsys):
        cfg = _valid_cfg()
        cfg["workflow"]["start_at"] = "phylo"
        cfg["workflow"]["stop_after"] = "hmmer"
        result = validate_config_only(cfg)
        assert result is False
        out = capsys.readouterr().out
        assert "start_at" in out or "stop_after" in out


# ── run_doctor ────────────────────────────────────────────────────────────────


class TestRunDoctor:
    def test_smoke_no_cfg(self, capsys):
        """run_doctor should not raise even when no tools are found."""
        run_doctor(cfg=None)
        out = capsys.readouterr().out
        assert "doctor" in out.lower() or "tool" in out.lower()

    def test_returns_bool(self):
        result = run_doctor(cfg=None)
        assert isinstance(result, bool)

    def test_with_valid_cfg(self, capsys):
        cfg = _valid_cfg()
        run_doctor(cfg=cfg)
        out = capsys.readouterr().out
        assert "Config validation" in out

    def test_all_required_tools_missing_returns_false(self, capsys):
        """When shutil.which always returns None, run_doctor returns False."""
        with patch("shutil.which", return_value=None):
            result = run_doctor(cfg=None)
        assert result is False

    def test_all_tools_present_returns_true(self, capsys):
        """When every tool is found, run_doctor returns True."""
        with patch("shutil.which", return_value="/usr/bin/fake"):
            result = run_doctor(cfg=None)
        assert result is True

    def test_python_packages_section(self, capsys):
        run_doctor(cfg=None)
        out = capsys.readouterr().out
        assert "Python packages" in out

    def test_required_tools_section(self, capsys):
        run_doctor(cfg=None)
        out = capsys.readouterr().out
        assert "Required tools" in out

    def test_optional_tools_section(self, capsys):
        run_doctor(cfg=None)
        out = capsys.readouterr().out
        assert "Optional tools" in out


# ── CLI integration (main.py) ─────────────────────────────────────────────────


class TestCLIIntegration:
    """Smoke-tests for the four new CLI flags via main.main()."""

    def _run_main(self, argv):
        """Run main() with patched sys.argv, capturing stdout and the exit code."""
        from phylofoundry.main import main

        with patch("sys.argv", ["phylofoundry"] + argv):
            with pytest.raises(SystemExit) as exc_info:
                main()
        return exc_info.value.code

    def test_list_steps_exits_zero(self, capsys):
        code = self._run_main(["--list-steps"])
        assert code == 0
        out = capsys.readouterr().out
        assert "prep" in out

    def test_plan_exits_zero_with_valid_args(self, capsys):
        code = self._run_main([
            "--plan",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code == 0
        out = capsys.readouterr().out
        assert "Execution plan" in out

    def test_validate_config_exits_zero_with_valid_args(self, capsys):
        code = self._run_main([
            "--validate-config",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code == 0

    def test_validate_config_exits_nonzero_when_missing_required(self, capsys):
        code = self._run_main(["--validate-config"])
        assert code != 0

    def test_doctor_exits_with_bool(self, capsys):
        code = self._run_main([
            "--doctor",
            "--faa_dir", "/tmp/faa",
            "--hmm_dir", "/tmp/hmm",
            "--outdir", "/tmp/out",
        ])
        assert code in (0, 1)
