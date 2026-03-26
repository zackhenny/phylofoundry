#!/usr/bin/env python3
"""migrate_json_to_yaml.py — Convert a PhyloFoundry JSON config to annotated YAML.

Usage
-----
    python scripts/migrate_json_to_yaml.py config.json
    python scripts/migrate_json_to_yaml.py config.json --output config/config.yaml

The script reads the existing JSON config, merges its values into the bundled
annotated YAML template (preserving all inline comments), and writes the result
to the specified output path (default: config/config.yaml alongside the input).

If no annotated template is found the script performs a plain JSON → YAML
conversion using ruamel.yaml's CommentedMap so that comments can be added
manually later.

Requirements
------------
    pip install ruamel.yaml
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _deep_update_commented(base, updates):
    """Recursively apply *updates* (plain dict) into *base* (CommentedMap).

    Nested dicts are merged rather than replaced so that inline comments on
    keys that already exist in *base* are preserved.
    """
    from ruamel.yaml.comments import CommentedMap, CommentedSeq

    for key, value in updates.items():
        if (
            isinstance(value, dict)
            and key in base
            and isinstance(base[key], (dict, CommentedMap))
        ):
            _deep_update_commented(base[key], value)
        elif isinstance(value, list) and key in base and isinstance(base[key], CommentedSeq):
            # Replace list contents while keeping any YAML-level comments on the key
            del base[key][:]
            base[key].extend(value)
        else:
            base[key] = value


def migrate(json_path: Path, output_path: Path, template_path: Path | None = None) -> None:
    """Perform the migration and write the output YAML file."""
    from ruamel.yaml import YAML
    from ruamel.yaml.comments import CommentedMap

    # ── Load the JSON source ────────────────────────────────────────────────
    print(f"[migrate] Reading JSON config: {json_path}")
    with json_path.open() as fh:
        json_cfg: dict = json.load(fh)

    yml = YAML()
    yml.preserve_quotes = True
    yml.default_flow_style = False
    yml.width = 120

    # ── Load the annotated template (preferred) ─────────────────────────────
    if template_path is None:
        # Try to find the bundled template relative to this script's location
        _script_dir = Path(__file__).parent
        _candidates = [
            _script_dir.parent / "config" / "config.yaml",
            _script_dir / ".." / "config" / "config.yaml",
        ]
        for c in _candidates:
            if c.resolve().exists():
                template_path = c.resolve()
                break

    if template_path and template_path.exists():
        print(f"[migrate] Using annotated template: {template_path}")
        with template_path.open() as fh:
            base = yml.load(fh)
        if base is None:
            base = CommentedMap()
        _deep_update_commented(base, json_cfg)
        header = (
            f"# Migrated from {json_path.name} by migrate_json_to_yaml.py\n"
            "# All values from the original JSON have been applied on top of the\n"
            "# annotated template.  Edit with ruamel.yaml to preserve comments.\n\n"
        )
    else:
        print(
            "[migrate] No annotated template found — performing plain JSON → YAML "
            "conversion.  Add comments manually or re-run after installing the "
            "bundled config/config.yaml template."
        )
        base = CommentedMap(json_cfg)
        header = (
            f"# Converted from {json_path.name} by migrate_json_to_yaml.py\n"
            "# Use ruamel.yaml to edit programmatically and preserve comments.\n\n"
        )

    # ── Write output ────────────────────────────────────────────────────────
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as fh:
        fh.write(header)
        yml.dump(base, fh)

    print(f"[migrate] Written YAML config: {output_path}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Convert a PhyloFoundry JSON config to annotated YAML.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "json_config",
        type=Path,
        help="Path to the existing config.json file to migrate.",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output path for the new YAML config "
            "(default: config/config.yaml relative to the JSON file's directory)."
        ),
    )
    ap.add_argument(
        "--template",
        type=Path,
        default=None,
        help=(
            "Path to the annotated YAML template to merge values into "
            "(default: auto-detected bundled config/config.yaml)."
        ),
    )
    args = ap.parse_args(argv)

    json_path: Path = args.json_config.resolve()
    if not json_path.exists():
        print(f"[migrate] ERROR: Input file not found: {json_path}", file=sys.stderr)
        sys.exit(1)

    if args.output is None:
        output_path = json_path.parent / "config" / "config.yaml"
    else:
        output_path = args.output.resolve()

    migrate(json_path, output_path, template_path=args.template)


if __name__ == "__main__":
    main()
