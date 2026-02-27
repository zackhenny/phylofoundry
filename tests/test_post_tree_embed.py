import re
from pathlib import Path

import pandas as pd

import phylofoundry.tasks.post as post_module


class _FakeNode:
    def __init__(self, name=None, length=0.0, children=None):
        self.name = name
        self.length = float(length)
        self.children = children or []
        self.parent = None
        self.id = None
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

    def distance(self, other):
        def ancestors(n):
            path = []
            cur = n
            d = 0.0
            while cur is not None:
                path.append((cur, d))
                d += cur.length
                cur = cur.parent
            return path

        a = ancestors(self)
        b = ancestors(other)
        da = {n: d for n, d in a}
        for n, db in b:
            if n in da:
                return da[n] + db
        return 0.0


class _FakeTree(_FakeNode):
    def assign_ids(self):
        for idx, n in enumerate(self.preorder(include_self=True), start=1):
            n.id = idx


def _synthetic_tree():
    left = _FakeNode("95", 1.0, [_FakeNode("A", 1.0), _FakeNode("B", 1.0)])
    right = _FakeNode("93", 1.0, [_FakeNode("C", 1.0), _FakeNode("D", 1.0)])
    root = _FakeTree("99", 0.0, [left, right])
    root.assign_ids()
    return root


def _write_fixture_tree_and_embeddings(root: Path, hmm: str = "HMMX"):
    tree_dir = root / "trees_iqtree"
    emb_dir = root / "embeddings"
    tree_dir.mkdir(parents=True, exist_ok=True)
    emb_dir.mkdir(parents=True, exist_ok=True)

    (tree_dir / f"{hmm}.treefile").write_text("((A:1,B:1)95:1,(C:1,D:1)93:1)99;\n")
    df = pd.DataFrame(
        [
            {"hmm": hmm, "tip": "A", "PC1": -5.2, "PC2": 0.1},
            {"hmm": hmm, "tip": "B", "PC1": -4.9, "PC2": -0.2},
            {"hmm": hmm, "tip": "C", "PC1": 5.1, "PC2": 0.1},
            {"hmm": hmm, "tip": "D", "PC1": 4.8, "PC2": -0.1},
        ]
    )
    df.to_csv(emb_dir / f"{hmm}.pca.tsv", sep="\t", index=False)
    return tree_dir, emb_dir


def test_detect_tree_embed_clades_synthetic_split(tmp_path, monkeypatch):
    tree_dir, emb_dir = _write_fixture_tree_and_embeddings(tmp_path)
    monkeypatch.setattr(post_module, "tree_load", lambda _: _synthetic_tree())

    post_cfg = {
        "embedtree_support_min": 80,
        "embedtree_min_size": 2,
        "embedtree_max_size": 100,
        "embedtree_top_k": 5,
        "embedtree_pcs": 2,
        "embedtree_distance": "euclidean",
        "embedtree_allow_nested": False,
    }
    clades, score_rows = post_module._detect_tree_embed_clades(str(tree_dir), str(emb_dir), post_cfg)

    assert clades
    assert score_rows
    assigned = {tip for tips in clades.values() for tip in tips}
    assert assigned.issubset({"A", "B", "C", "D"})


def test_run_post_tree_embed_writes_detected_clades(tmp_path, monkeypatch):
    results_dir = tmp_path / "results"
    summary_dir = results_dir / "summary"
    post_dir = summary_dir / "post_scikitbio"
    clipkit_dir = results_dir / "alignments_clipkit"
    aln_dir = results_dir / "alignments"
    tree_dir, _ = _write_fixture_tree_and_embeddings(results_dir, hmm="HMMY")

    summary_dir.mkdir(parents=True, exist_ok=True)
    post_dir.mkdir(parents=True, exist_ok=True)
    clipkit_dir.mkdir(parents=True, exist_ok=True)
    aln_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(post_module, "tree_load", lambda _: _synthetic_tree())

    cfg = {
        "inputs": {"gtdb_dir": None, "taxonomy_file": None},
        "post": {
            "detect_clades_method": "tree_embed",
            "embedtree_support_min": 0,
            "embedtree_min_size": 2,
            "embedtree_top_k": 3,
            "embedtree_pcs": 2,
            "embedtree_distance": "euclidean",
            "compute_conservation": False,
            "compute_kl": False,
        },
    }

    post_module.run_post(
        cfg,
        str(tree_dir),
        str(clipkit_dir),
        str(aln_dir),
        str(post_dir),
        str(summary_dir),
        hmm_keep=None,
        force=True,
    )

    detected_fp = summary_dir / "detected_clades.tsv"
    scores_fp = summary_dir / "node_scores.embedtree.tsv"
    assert detected_fp.exists()
    assert scores_fp.exists()

    detected = pd.read_csv(detected_fp, sep="\t")
    assert {"clade_name", "tip"}.issubset(detected.columns)
    assert detected["clade_name"].map(lambda x: bool(re.match(r"embed__C\d+", str(x)))).all()
    assert set(detected["tip"]).issubset({"A", "B", "C", "D"})
