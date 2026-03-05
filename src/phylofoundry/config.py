from __future__ import annotations

import argparse
import json
import os
import shutil
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError

from .constants import DEFAULT_CONFIG, STEPS
from .utils.helpers import load_json_config


class InputsConfig(BaseModel):
    faa_dir: str | None = Field(default=None, description="Directory or single .faa file")
    hmm_input: str | None = Field(default=None, description="Directory or single .hmm file")
    cds_dir: str | None = None


class OutputConfig(BaseModel):
    outdir: str | None = Field(default=None, description="Pipeline output directory")


class ResourceConfig(BaseModel):
    cpu: int = Field(default=8, ge=1, description="CPU threads")


class WorkflowConfig(BaseModel):
    start_at: str | None = None
    stop_after: str | None = None
    force: bool = False
    hmm_manifest: str | None = None
    preset: str | None = Field(default=None, description="Workflow preset: phylo_only|plm_hypothesis|plm_plus_selection")
    max_failure_rate: float = Field(default=0.2, ge=0.0, le=1.0)


class EmbeddingsConfig(BaseModel):
    enabled: bool = False
    backend: str = "esm"
    model: str = "esm2_t33_650M_UR50D"
    device: str = "cuda"
    batch_size: int = Field(default=8, ge=1)
    write_full_vectors: bool = False


class RegimeShiftConfig(BaseModel):
    enable: bool = False
    metric: str = Field(default="centroid", description="centroid|energy|mmd")
    min_support: float = 0.0
    min_size: int = Field(default=3, ge=2)
    n_permutations: int = Field(default=200, ge=0)
    alpha: float = Field(default=0.05, gt=0.0, le=1.0)
    require_monophyly: bool = True


class HAConfig(BaseModel):
    enabled: bool = False
    mode: str = Field(default="loc", description="loc|middle")
    pooling_used: str = Field(default="mean", description="mean|max|both")
    loc_theta_target_deg: float = 90.0
    loc_break_adjust: int = -1
    layers: list[int] | None = None
    call_mode: str = Field(default="loc_break", description="loc_break|percentile|topk")
    percentile: float = Field(default=0.95, gt=0.0, lt=1.0)
    topk: int = Field(default=20, ge=1)
    max_logged_failures: int = Field(default=5, ge=1)


class EvidenceJoinConfig(BaseModel):
    enable: bool = False
    classification_thresholds: dict[str, float] = Field(default_factory=lambda: {"delta_ha": 0.2, "js": 0.15, "meme_p": 0.1})
    grading_thresholds: dict[str, float] = Field(default_factory=lambda: {"A": 0.8, "B": 0.5})


class QCConfig(BaseModel):
    enable: bool = True
    per_hmm: bool = True
    combined: bool = True
    max_hmms_to_plot: int = Field(default=100, ge=1)
    dpi: int = Field(default=150, ge=72)


class RootConfig(BaseModel):
    model_config = ConfigDict(extra="allow")

    inputs: InputsConfig = Field(default_factory=InputsConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    resources: ResourceConfig = Field(default_factory=ResourceConfig)
    workflow: WorkflowConfig = Field(default_factory=WorkflowConfig)
    embeddings: EmbeddingsConfig = Field(default_factory=EmbeddingsConfig)
    regime_shift: RegimeShiftConfig = Field(default_factory=RegimeShiftConfig)
    ha: HAConfig = Field(default_factory=HAConfig)
    evidence_join: EvidenceJoinConfig = Field(default_factory=EvidenceJoinConfig)
    qc: QCConfig = Field(default_factory=QCConfig)


def deep_update(base: dict, updates: dict) -> dict:
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base


def load_config_file(path: str) -> dict:
    p = Path(path)
    if p.suffix.lower() in {".yaml", ".yml"}:
        with open(p) as f:
            return yaml.safe_load(f) or {}
    return load_json_config(path)


def apply_set_overrides(cfg: dict, overrides: list[str] | None):
    if not overrides:
        return
    for raw in overrides:
        if "=" not in raw:
            raise SystemExit(f"Invalid --set override '{raw}'. Expected section.key=value")
        key, val = raw.split("=", 1)
        tgt = cfg
        keys = key.split(".")
        for part in keys[:-1]:
            tgt = tgt.setdefault(part, {})
        try:
            tgt[keys[-1]] = json.loads(val)
        except json.JSONDecodeError:
            tgt[keys[-1]] = val


def apply_workflow_preset(cfg: dict):
    preset = cfg.get("workflow", {}).get("preset")
    if not preset:
        return
    if preset == "phylo_only":
        cfg.setdefault("embeddings", {})["enabled"] = False
        cfg.setdefault("regime_shift", {})["enable"] = False
        cfg.setdefault("ha", {})["enabled"] = False
        cfg.setdefault("hyphy", {})["enabled"] = False
    elif preset == "plm_hypothesis":
        cfg.setdefault("embeddings", {})["enabled"] = True
        cfg.setdefault("regime_shift", {})["enable"] = True
        cfg.setdefault("ha", {})["enabled"] = True
        cfg.setdefault("evidence_join", {})["enable"] = True
    elif preset == "plm_plus_selection":
        cfg.setdefault("embeddings", {})["enabled"] = True
        cfg.setdefault("regime_shift", {})["enable"] = True
        cfg.setdefault("ha", {})["enabled"] = True
        cfg.setdefault("codon", {})["enabled"] = True
        cfg.setdefault("hyphy", {})["enabled"] = True
        cfg.setdefault("evidence_join", {})["enable"] = True
    else:
        raise SystemExit(f"Unknown workflow.preset '{preset}'")


def resolve_config(args: argparse.Namespace) -> dict:
    if getattr(args, "dump_default_config", False):
        print(json.dumps(DEFAULT_CONFIG, indent=2, sort_keys=True))
        return None

    cfg = deepcopy(DEFAULT_CONFIG)
    if getattr(args, "config", None):
        deep_update(cfg, load_config_file(args.config))

    for src, sec, key in [
        ("faa_dir", "inputs", "faa_dir"),
        ("hmm_dir", "inputs", "hmm_input"),
        ("outdir", "output", "outdir"),
    ]:
        v = getattr(args, src, None)
        if v is not None:
            cfg[sec][key] = v
    if getattr(args, "cpu", None) is not None:
        cfg["resources"]["cpu"] = int(args.cpu)
    if getattr(args, "cpu", None) is None and (os.environ.get("SLURM_CPUS_PER_TASK") or "").isdigit():
        cfg["resources"]["cpu"] = int(os.environ["SLURM_CPUS_PER_TASK"])

    if getattr(args, "start_at", None) is not None:
        cfg["workflow"]["start_at"] = args.start_at
    if getattr(args, "stop_after", None) is not None:
        cfg["workflow"]["stop_after"] = args.stop_after
    if getattr(args, "force", False):
        cfg["workflow"]["force"] = True

    apply_set_overrides(cfg, getattr(args, "set_values", None))
    apply_workflow_preset(cfg)

    current_bin = cfg.get("phylo", {}).get("iqtree_bin", "iqtree")
    if current_bin == "iqtree" and not shutil.which("iqtree"):
        for cand in ["iqtree2", "iqtree3"]:
            if shutil.which(cand):
                cfg.setdefault("phylo", {})["iqtree_bin"] = cand
                break
    return cfg


def validate_config(cfg: dict) -> dict:
    try:
        model = RootConfig.model_validate(cfg)
    except ValidationError as e:
        raise SystemExit(f"Config validation failed:\n{e}") from e

    data = model.model_dump()
    if not data["inputs"]["faa_dir"] or not data["inputs"]["hmm_input"] or not data["output"]["outdir"]:
        raise SystemExit("Config must specify inputs.faa_dir, inputs.hmm_input, output.outdir (or pass via CLI).")
    for k in ["start_at", "stop_after"]:
        v = data["workflow"].get(k)
        if v and v not in STEPS:
            raise SystemExit(f"workflow.{k} must be one of {STEPS}")
    return data


def config_template(mode: str = "minimal") -> str:
    if mode == "full":
        return yaml.safe_dump(DEFAULT_CONFIG, sort_keys=False)
    minimal = {
        "inputs": {"faa_dir": "path/to/faa", "hmm_input": "path/to/hmms"},
        "output": {"outdir": "phylofoundry_out"},
        "workflow": {"preset": "plm_hypothesis"},
        "resources": {"cpu": 8},
    }
    return yaml.safe_dump(minimal, sort_keys=False)


def config_explain(path: str | None = None) -> str:
    sections = {
        "inputs": "Input data sources.",
        "output": "Output location.",
        "workflow": "Execution controls + presets.",
        "embeddings": "PLM embedding generation.",
        "regime_shift": "Embedding-based clade/branch shift detection.",
        "ha": "High-attention residue calling.",
        "evidence_join": "Integrative hypothesis classification.",
        "qc": "QC plot generation controls.",
    }
    if path:
        cfg = load_config_file(path)
        return yaml.safe_dump(cfg, sort_keys=False)
    return "\n".join(f"{k}: {v}" for k, v in sections.items())
