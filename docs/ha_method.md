# HA-site method (independent module)

PhyloFoundry now runs HA-site analysis as an independent stage (`phylofoundry ha`) aligned with Helix-Research-Lab/ESM_highAttention.

## Reference-aligned algorithm

1. Extract `results["attentions"]` from ESM output (`[B, layers, heads, T, T]`).
2. Pool over heads with both `mean` and `max`.
3. Stack to `combined_attn_tensor` with shape `[2, B, layers, T, T]` and save via `torch.save`.
4. For HA inference use mean channel by default: `data[0,0,layer,1:-1,1:-1]`.
5. Compute received attention by column sum: `sum(attn_matrix, dim=0)`.
6. Normalize each layer vector by its max.
7. LoC mode:
   - sort descending,
   - fit 2-piece linear model (`pwlf`) on `x = linspace(0,1,L)`,
   - take breakpoint `breaks[1]` as HA fraction,
   - `n_ha = floor(break_frac * L) + loc_break_adjust`.
8. Compute `theta` from piecewise slopes and pick LoC layer closest to `loc_theta_target_deg` (default 90°).

## Key parameters mapping

- `ha.mode`: `loc|middle`
- `ha.pooling_used`: `mean|max|both`
- `ha.call_mode`: `loc_break|percentile|topk`
- `ha.loc_theta_target_deg`: LoC target angle
- `ha.loc_break_adjust`: additive correction to breakpoint-derived site count
- `ha.layers`: optional explicit layers

## Outputs

Per HMM:
- `summary/ha_sites.tsv`
- `summary/ha_counts.tsv`
- `summary/loc_layers.tsv`
- `ha/heatmaps/<seq_id>.npy`
- `ha/attn/<seq_id>.pt`
- `qc/<hmm_id>/ha/*.png`

## Interpretation caveat

HA sites are model-salient positions from the PLM attention mechanism; they should be interpreted as hypotheses and combined with evolutionary/structural evidence.
