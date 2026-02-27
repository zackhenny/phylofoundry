from pathlib import Path

import pandas as pd

import phylofoundry.tasks.hyphy as hyphy


class FakeNode:
    def __init__(self, name=None, children=None):
        self.name = name
        self.children = children or []
        self.parent = None
        for c in self.children:
            c.parent = self

    def is_tip(self):
        return len(self.children) == 0

    def tips(self):
        if self.is_tip():
            return [self]
        out = []
        for c in self.children:
            out.extend(c.tips())
        return out

    def preorder(self, include_self=True):
        nodes = [self] if include_self else []
        for c in self.children:
            nodes.extend(c.preorder(include_self=True))
        return nodes

    def copy(self):
        return _clone(self)

    def write(self, out_fp, format="newick"):
        assert format == "newick"

        def emit(node):
            if node.is_tip():
                return node.name
            kids = ",".join(emit(c) for c in node.children)
            label = node.name or ""
            return f"({kids}){label}"

        Path(out_fp).write_text(emit(self) + ";\n")


class FakeTree(FakeNode):
    def lca(self, names):
        names = set(names)

        def recurse(node):
            tipset = {t.name for t in node.tips()}
            if names.issubset(tipset):
                for child in node.children:
                    cset = {t.name for t in child.tips()}
                    if names.issubset(cset):
                        return recurse(child)
                return node
            return None

        return recurse(self)


def _clone(node):
    copied = FakeTree(node.name) if isinstance(node, FakeTree) else FakeNode(node.name)
    copied.children = [_clone(c) for c in node.children]
    for c in copied.children:
        c.parent = copied
    return copied


def _synthetic_tree():
    left_fg = FakeNode(None, [FakeNode("t1"), FakeNode("t2"), FakeNode("t3")])
    left_bg = FakeNode(None, [FakeNode("t4")])
    left = FakeNode(None, [left_fg, left_bg])
    right = FakeNode(None, [FakeNode("t5"), FakeNode("t6")])
    return FakeTree(None, [left, right])


def test_build_labeled_tree_labels_and_tip_names():
    tree = _synthetic_tree()
    relax_tree = hyphy._build_labeled_tree_for_clade(
        tree,
        ["t1", "t2", "t3"],
        test_name="RELAX",
        label_mode="stem",
        relax_label_reference=True,
    )
    tip_names = {n.name for n in relax_tree.tips()}
    assert tip_names == {"t1", "t2", "t3", "t4", "t5", "t6"}

    internal_labels = [n.name or "" for n in relax_tree.preorder() if not n.is_tip()]
    assert any("{test}" in x for x in internal_labels)
    assert any("{reference}" in x for x in internal_labels)

    fg_tree = hyphy._build_labeled_tree_for_clade(
        tree,
        ["t1", "t2", "t3"],
        test_name="aBSREL",
        label_mode="crown",
        relax_label_reference=True,
    )
    fg_labels = [n.name or "" for n in fg_tree.preorder() if not n.is_tip()]
    assert any("{FG}" in x for x in fg_labels)


def test_run_hyphy_clade_aware_builds_commands_and_outputs(tmp_path, monkeypatch):
    results_dir = tmp_path / "results"
    codon_dir = results_dir / "codon_alignments"
    tree_dir = results_dir / "trees_iqtree"
    summary_dir = results_dir / "summary"
    hyphy_dir = summary_dir / "hyphy"
    codon_dir.mkdir(parents=True)
    tree_dir.mkdir(parents=True)
    hyphy_dir.mkdir(parents=True)

    (codon_dir / "HMM1.codon.fasta").write_text(">t1\nATGATG\n>t2\nATGATG\n>t3\nATGATG\n>t4\nATGATG\n>t5\nATGATG\n>t6\nATGATG\n")
    (tree_dir / "HMM1.treefile").write_text("dummy;\n")

    pd.DataFrame(
        [
            {"clade_name": "cladeA", "tip": "t1"},
            {"clade_name": "cladeA", "tip": "t2"},
            {"clade_name": "cladeA", "tip": "t3"},
            {"clade_name": "cladeA", "tip": "missing"},
        ]
    ).to_csv(summary_dir / "detected_clades.tsv", sep="\t", index=False)

    monkeypatch.setattr(hyphy, "_import_post_helpers", lambda: (lambda fp: {"cladeA": ["t1", "t2", "t3", "missing"]}, lambda fp: _synthetic_tree()))

    calls = []

    def fake_wrapper(hyphy_bin, test_name, codon_aln_fp, tree_fp, out_json, extra_args=None, log_fp=None):
        calls.append((test_name, tree_fp, out_json, extra_args, log_fp))
        Path(out_json).parent.mkdir(parents=True, exist_ok=True)
        Path(out_json).write_text("{}")
        if log_fp:
            Path(log_fp).parent.mkdir(parents=True, exist_ok=True)
            Path(log_fp).write_text("ok")
        return True, [hyphy_bin, test_name]

    monkeypatch.setattr(hyphy, "run_hyphy_wrapper", fake_wrapper)

    cfg = {
        "hyphy": {
            "hyphy_bin": "hyphy",
            "hyphy_tests": "RELAX,aBSREL,BUSTED,MEME",
            "use_detected_clades": True,
            "min_clade_size": 3,
            "label_mode": "crown",
            "relax_label_reference": True,
            "hyphy_args": {
                "MEME": ["--branches", "All"],
                "aBSREL": ["--branches", "All"],
                "BUSTED": ["--branches", "All"],
                "RELAX": ["--test", "Test", "--reference", "Reference"],
            },
        }
    }

    hyphy.run_hyphy(cfg, str(codon_dir), str(tree_dir), str(hyphy_dir), hmm_keep=None, force=True)

    assert any(c[0].upper() == "RELAX" for c in calls)
    assert any(c[0].upper() == "ABSREL" for c in calls)
    assert any(c[0].upper() == "BUSTED" for c in calls)

    relax_call = next(c for c in calls if c[0].upper() == "RELAX")
    assert "--test" in relax_call[3] and "test" in relax_call[3]
    assert "--reference" in relax_call[3] and "reference" in relax_call[3]

    absrel_call = next(c for c in calls if c[0].upper() == "ABSREL")
    assert "--branches" in absrel_call[3] and "FG" in absrel_call[3]

    assert Path(results_dir / "trees_labeled" / "HMM1" / "cladeA.RELAX.nwk").exists()
    assert Path(results_dir / "trees_labeled" / "HMM1" / "cladeA.ABSREL.nwk").exists()
    assert Path(hyphy_dir / "clade_runs.tsv").exists()
