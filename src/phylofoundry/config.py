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
    faa_dir: str | None = None
    hmm_input: str | None = None
    cds_dir: str | None = None


class OutputConfig(BaseModel):
    outdir: str | None = None


class ResourceConfig(BaseModel):
    cpu: int = Field(default=8, ge=1)


class WorkflowConfig(BaseModel):
    start_at: str | None = None
    stop_after: str | None = None
    force: bool = False
    hmm_manifest: str | None = None


class EmbeddingsConfig(BaseModel):
    enabled: bool = False
    device: str = "cuda"
    batch_size: int = Field(default=8, ge=1)
    write_full_vectors: bool = False


class DiscoverConfig(BaseModel):
    enabled: bool = False
    kmer_size: int = Field(default=5, ge=3)
    top_n_peaks: int = Field(default=20, ge=1)
    attention_layers: int = Field(default=4, ge=1)


class PhyloConfig(BaseModel):
    combined_tree: bool = False
    iqtree_bin: str = "iqtree"


class RootConfig(BaseModel):
    model_config = ConfigDict(extra="allow")

    inputs: InputsConfig = Field(default_factory=InputsConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    resources: ResourceConfig = Field(default_factory=ResourceConfig)
    workflow: WorkflowConfig = Field(default_factory=WorkflowConfig)
    embeddings: EmbeddingsConfig = Field(default_factory=EmbeddingsConfig)
    discover: DiscoverConfig = Field(default_factory=DiscoverConfig)
    phylo: PhyloConfig = Field(default_factory=PhyloConfig)


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
            data = yaml.safe_load(f) or {}
        return data
    return load_json_config(path)


def apply_set_overrides(cfg: dict, overrides: list[str] | None):
    if not overrides:
        return
    for raw in overrides:
        if "=" not in raw:
            raise SystemExit(f"Invalid --set override '{raw}'. Expected section.key=value")
        key, val = raw.split("=", 1)
        keys = key.split(".")
        tgt = cfg
        for part in keys[:-1]:
            tgt = tgt.setdefault(part, {})
        try:
            parsed = json.loads(val)
        except json.JSONDecodeError:
            parsed = val
        tgt[keys[-1]] = parsed


def resolve_config(args: argparse.Namespace) -> dict:
    if getattr(args, "dump_default_config", False):
        print(json.dumps(DEFAULT_CONFIG, indent=2, sort_keys=True))
        return None

    cfg = deepcopy(DEFAULT_CONFIG)
    if getattr(args, "config", None):
        deep_update(cfg, load_config_file(args.config))

    if getattr(args, "faa_dir", None) is not None:
        cfg["inputs"]["faa_dir"] = args.faa_dir
    if getattr(args, "hmm_dir", None) is not None:
        cfg["inputs"]["hmm_input"] = args.hmm_dir
    if getattr(args, "outdir", None) is not None:
        cfg["output"]["outdir"] = args.outdir
    if getattr(args, "cpu", None) is not None:
        cfg["resources"]["cpu"] = int(args.cpu)

    if getattr(args, "cpu", None) is None:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        if slurm_cpus and str(slurm_cpus).isdigit():
            cfg["resources"]["cpu"] = int(slurm_cpus)

    if getattr(args, "start_at", None) is not None:
        cfg["workflow"]["start_at"] = args.start_at
    if getattr(args, "stop_after", None) is not None:
        cfg["workflow"]["stop_after"] = args.stop_after
    if getattr(args, "force", False):
        cfg["workflow"]["force"] = True

    apply_set_overrides(cfg, getattr(args, "set_values", None))

    current_bin = cfg["phylo"].get("iqtree_bin", "iqtree")
    if current_bin == "iqtree" and not shutil.which("iqtree"):
        for cand in ["iqtree2", "iqtree3"]:
            if shutil.which(cand):
                cfg["phylo"]["iqtree_bin"] = cand
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

    start_at = data["workflow"].get("start_at")
    stop_after = data["workflow"].get("stop_after")
    if start_at and start_at not in STEPS:
        raise SystemExit(f"workflow.start_at must be one of {STEPS}")
    if stop_after and stop_after not in STEPS:
        raise SystemExit(f"workflow.stop_after must be one of {STEPS}")
    return data


def config_template() -> str:
    minimal = {
        "inputs": {"faa_dir": "path/to/faa", "hmm_input": "path/to/hmms"},
        "output": {"outdir": "phylofoundry_out"},
        "resources": {"cpu": 8},
        "embeddings": {"enabled": False, "device": "cpu"},
    }
    return yaml.safe_dump(minimal, sort_keys=False)


def config_explain() -> str:
    return (
        "inputs: input data sources (FAA + HMM are required)\n"
        "output: output directory\n"
        "resources: cpu / runtime resources\n"
        "workflow: start_at/stop_after for partial runs\n"
        "embeddings: ESM/transformer embeddings behavior\n"
        "discover: comparative attention-based motif discovery\n"
    )
