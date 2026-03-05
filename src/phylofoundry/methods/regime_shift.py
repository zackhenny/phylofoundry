from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class Node:
    name: str | None = None
    children: list["Node"] = field(default_factory=list)
    support: float | None = None
    node_id: str | None = None

    def is_tip(self) -> bool:
        return len(self.children) == 0


def _parse_newick_simple(newick: str) -> Node:
    stack: list[Node] = []
    cur = Node()
    token = ""
    for ch in newick.strip().rstrip(";"):
        if ch == "(":
            stack.append(cur)
            n = Node()
            cur.children.append(n)
            cur = n
        elif ch == ",":
            if token.strip():
                cur.name = token.split(":")[0].strip()
                token = ""
            parent = stack[-1]
            sib = Node()
            parent.children.append(sib)
            cur = sib
        elif ch == ")":
            if token.strip():
                cur.name = token.split(":")[0].strip()
                token = ""
            cur = stack.pop()
        else:
            token += ch
    return cur


def _tip_names(n: Node) -> list[str]:
    if n.is_tip():
        return [n.name] if n.name else []
    out = []
    for c in n.children:
        out.extend(_tip_names(c))
    return out


def _iter_internal(n: Node):
    if not n.is_tip():
        yield n
    for c in n.children:
        yield from _iter_internal(c)


def _centroid(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a.mean(0) - b.mean(0)))


def _energy(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.spatial.distance import cdist

    return float(2 * cdist(a, b).mean() - cdist(a, a).mean() - cdist(b, b).mean())


def _mmd(a: np.ndarray, b: np.ndarray) -> float:
    gamma = 1.0 / max(a.shape[1], 1)

    def k(x, y):
        d = ((x[:, None, :] - y[None, :, :]) ** 2).sum(-1)
        return np.exp(-gamma * d)

    return float(k(a, a).mean() + k(b, b).mean() - 2 * k(a, b).mean())


def _score(a: np.ndarray, b: np.ndarray, metric: str) -> float:
    if metric == "energy":
        return _energy(a, b)
    if metric == "mmd":
        return _mmd(a, b)
    return _centroid(a, b)


def run_regime_shift(cfg: dict, tree_dir: str, emb_dir: str, summary_dir: str, hmm_keep=None):
    rs = cfg.get("regime_shift", {})
    metric = rs.get("metric", "centroid")
    min_size = int(rs.get("min_size", 3))
    n_perm = int(rs.get("n_permutations", 200))

    rows = []
    ann = []
    for pca_fp in sorted(glob.glob(os.path.join(emb_dir, "*.pca.tsv"))):
        hmm = os.path.basename(pca_fp).replace(".pca.tsv", "")
        if hmm_keep and hmm not in hmm_keep:
            continue
        tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
        if not os.path.exists(tree_fp):
            raise SystemExit(f"RegimeShift missing required tree: {tree_fp}")
        df = pd.read_csv(pca_fp, sep="\t")
        feats = [c for c in df.columns if c.startswith("PC")]
        emap = {r["tip"]: r[feats].to_numpy(dtype=float) for _, r in df.iterrows()}

        root = _parse_newick_simple(open(tree_fp).read())
        node_i = 0
        for n in _iter_internal(root):
            if len(n.children) < 2:
                continue
            left = [t for t in _tip_names(n.children[0]) if t in emap]
            right = [t for t in _tip_names(n.children[1]) if t in emap]
            if len(left) < min_size or len(right) < min_size:
                continue
            A = np.vstack([emap[t] for t in left])
            B = np.vstack([emap[t] for t in right])
            obs = _score(A, B, metric)
            p = np.nan
            if n_perm > 0:
                pooled = np.vstack([A, B])
                n_a = len(A)
                perm = []
                rng = np.random.default_rng(0)
                for _ in range(n_perm):
                    idx = rng.permutation(len(pooled))
                    perm.append(_score(pooled[idx[:n_a]], pooled[idx[n_a:]], metric))
                p = float((np.sum(np.asarray(perm) >= obs) + 1) / (n_perm + 1))
            node_i += 1
            nid = f"{hmm}_N{node_i}"
            rows.append({"hmm_id": hmm, "node_id": nid, "support": np.nan, "n_left": len(left), "n_right": len(right), "shift_score": obs, "p_perm": p, "embedding_distance_metric": metric, "notes": ""})
            ann.append({"hmm_id": hmm, "node_id": nid, "shift_score": obs})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(summary_dir, "regime_shifts.tsv"), sep="\t", index=False)
    nwk_fp = os.path.join(summary_dir, "regime_shift_nodes.nwk")
    with open(nwk_fp, "w") as f:
        for r in ann:
            f.write(f"[{r['hmm_id']}|{r['node_id']}|shift={r['shift_score']:.6f}]\n")
    return out
