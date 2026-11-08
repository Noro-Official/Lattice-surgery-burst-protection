"""Continuous-history d=13 memory trajectories with round-resolved burst noise (2026-10-01; stage B pilot).

One trajectory = `prefix` quiet rounds (background only), then an unannounced burst at round `prefix`, then recovery.
Round r >= prefix sits at cycle (r - prefix) after the impact; its channel-B site rates and its state-conditioned
readout rates (eps0, eps1) are evaluated at that cycle (no bin averaging).  Same template, slots, weights and
readout treatment as the block campaign (`pacer_readout_model.fill_timed_readout` + `sample_conditioned`), so a
trajectory with constant per-round values is IDENTICAL to a block circuit (checked by `selfcheck`).
"""
from __future__ import annotations

from typing import Any

import numpy as np
import stim

from . import pacer_readout_model as rm
from . import pacer_timed_channel as tc


def split_rounds(template: Any):
    """(chunks, groups, round_of_group): tc.split_timed plus the syndrome round each slot group belongs to.

    A round ends at its MR; the post-reset slot directly after an MR belongs to the round just measured; the final
    data-measurement slots belong to round R (one past the last syndrome round)."""
    chunks, groups, rid, current, group = [], [], [], stim.Circuit(), []
    r, after_mr = 0, False
    grp_round = None
    for inst in template:
        if inst.name == "PAULI_CHANNEL_1" and getattr(inst, "tag", ""):
            if not group:
                chunks.append(current)
                current = stim.Circuit()
                grp_round = (r - 1) if (after_mr and inst.tag == "anc_post_reset") else r
            for t in inst.targets_copy():
                group.append((inst.tag, int(t.value)))
        else:
            if group:
                groups.append(group); rid.append(grp_round); group = []
            current.append(inst)
            if inst.name == "MR":
                r += 1
                after_mr = True
            elif inst.name not in ("TICK", "DETECTOR", "QUBIT_COORDS", "SHIFT_COORDS", "X_ERROR", "DEPOLARIZE1"):
                after_mr = False
    if group:
        groups.append(group); rid.append(grp_round)
    chunks.append(current)
    return chunks, groups, rid


def fill_rounds(chunks, groups, rid, layout, g_data, q_data, g_anc, q_anc, p_reset, kind_weight: dict | None = None,
                p_meas=None):
    """Round-resolved fill.  Arrays are [R + 1, n] (row R: the final data measurement); p_reset is [R + 1, n_anc].
    The anc_meas slot carries no burst readout excess (state-conditioned flips are applied to the records)."""
    kw = kind_weight or rm.physical_kind_weight()
    dpos = {int(q): i for i, q in enumerate(layout.data_ids)}
    apos = {int(q): i for i, q in enumerate(layout.ancilla_ids)}
    out = stim.Circuit()
    nrow = len(g_data) - 1
    for i, chunk in enumerate(chunks):
        out += chunk
        if i < len(groups):
            r = min(rid[i], nrow)
            parts = []
            for kind, qb in groups[i]:
                if kind == "anc_meas":
                    # p_meas given (decoder-model circuits only): symmetric readout flip X_ERROR((eps0 + eps1) / 2)
                    if p_meas is not None and float(p_meas[r, apos[qb]]) > 0:
                        parts.append(f"X_ERROR({float(p_meas[r, apos[qb]])!r}) {qb}")
                    continue
                if kind == "anc_post_reset":
                    p = float(p_reset[r, apos[qb]])
                    if p > 0:
                        parts.append(f"X_ERROR({p!r}) {qb}")
                    continue
                if qb in dpos:
                    g, q = float(g_data[r, dpos[qb]]), float(q_data[r, dpos[qb]])
                else:
                    g, q = float(g_anc[r, apos[qb]]), float(q_anc[r, apos[qb]])
                if g > 0 or q > 0:
                    pp = tc.layer_pauli(g, q, kw[kind])
                    if any(pp):
                        parts.append("PAULI_CHANNEL_1({!r},{!r},{!r}) {}".format(pp[0], pp[1], pp[2], qb))
            if parts:
                out += stim.Circuit("\n".join(parts))
    return out


def trajectory_channels(x0_event, design, pstar, layout, adapter, rounds: int, prefix: int,
                        decoupling: str = "standard"):
    """Per-round arrays for one event/design: (g_data, q_data, g_anc, q_anc, p_reset, eps0, eps1), rows 0..R.
    Rows < prefix (and nothing else) are burst-free.  Row R (final data measurement) uses cycle R - prefix.
    decoupling: bridge Platform decoupling ("standard" = Willow circuit (i), the historical default; "echo" = circuit (iii))."""
    from ..physics.pacer_qp_bridge_v1 import channels_from_site_design
    nd, na = len(layout.data_ids), len(layout.ancilla_ids)
    R = rounds
    cyc = np.arange(prefix, R + 1) - prefix
    times = (cyc + 0.5) * pstar.cycle_time_s
    ch = channels_from_site_design(np.asarray(x0_event, float), times, pstar, gap_difference_hz=np.asarray(design["gaps"], float),
                                   trapping_hz=np.asarray(design["traps"], float), decoupling=decoupling)
    gd = np.zeros((R + 1, nd)); qd = np.zeros((R + 1, nd)); ga = np.zeros((R + 1, na)); qa = np.zeros((R + 1, na))
    pr = np.zeros((R + 1, na)); e0 = np.zeros((R + 1, na)); e1 = np.zeros((R + 1, na))
    pdd, pda = np.asarray(ch.pauli_data), np.asarray(ch.pauli_ancilla)
    gdn, gup = rm.rates_from_event_channels(ch, slice(0, len(cyc)), nd, pstar.cycle_time_s)
    probs = adapter.probabilities(gdn, gup)
    for k in range(len(cyc)):
        r = prefix + k
        gd[r], qd[r] = tc.invert(pdd[k, :nd, :])
        ga[r], qa[r] = tc.invert(pda[k, nd:, :])
    e0[prefix:], e1[prefix:] = probs["flip_excess"][..., 0], probs["flip_excess"][..., 1]
    pr[prefix:] = probs["reset_excess"].mean(-1)
    return gd, qd, ga, qa, pr, e0, e1


def sample_trajectory(circuit, shots: int, seed: int, rec_idx, anc_pos, eps0, eps1, rounds: int,
                      return_physical: bool = False):
    """Raw records + state-conditioned flips with per-round rates.  eps0/eps1: [R + 1, n_anc] (rows 0..R-1 used).
    Returns (det, obs, y) with y the REPORTED ancilla records [shots, R, n_anc] (what a decoder may see)."""
    meas = circuit.compile_sampler(seed=seed).sample(int(shots))
    m = meas[:, rec_idx]
    n_anc = len(rec_idx) // rounds
    rr = np.repeat(np.arange(rounds), n_anc)                    # records are round-major over ancillas
    e0, e1 = np.asarray(eps0)[rr, anc_pos], np.asarray(eps1)[rr, anc_pos]
    rng = np.random.default_rng([int(seed) % (2 ** 63), 1])
    flip = rng.random(m.shape) < np.where(m, e1, e0)
    y = m ^ flip
    meas[:, rec_idx] = y
    det, obs = circuit.compile_m2d_converter().convert(measurements=meas, separate_observables=True)
    if return_physical:      # m = the physical ancilla bit immediately BEFORE the burst readout error (evaluation only)
        return det, obs, y.reshape(len(y), rounds, n_anc), m.reshape(len(m), rounds, n_anc)
    return det, obs, y.reshape(len(y), rounds, n_anc)


# --- fast path for NEW cohorts (2026-10-02): same circuits and the same sampling distribution, less work ---------
def fill_rounds_pair(chunks, groups, rid, layout, g_data, q_data, g_anc, q_anc, p_reset, p_meas=None, kind_weight: dict | None = None):
    """One pass over the template producing (trajectory circuit, decoder-model circuit or None).
    The model circuit is the trajectory circuit plus X_ERROR(p_meas) in every anc_meas slot (state-averaged readout).
    Identical, instruction for instruction, to fill_rounds(...) and fill_rounds(..., p_meas=p_meas)."""
    kw = kind_weight or rm.physical_kind_weight()
    dpos = {int(q): i for i, q in enumerate(layout.data_ids)}
    apos = {int(q): i for i, q in enumerate(layout.ancilla_ids)}
    out = stim.Circuit()
    mod = stim.Circuit() if p_meas is not None else None
    nrow = len(g_data) - 1
    for i, chunk in enumerate(chunks):
        out += chunk
        if mod is not None:
            mod += chunk
        if i < len(groups):
            r = min(rid[i], nrow)
            parts, mparts = [], []
            for kind, qb in groups[i]:
                if kind == "anc_meas":
                    if mod is not None and float(p_meas[r, apos[qb]]) > 0:
                        mparts.append(f"X_ERROR({float(p_meas[r, apos[qb]])!r}) {qb}")
                    continue
                if kind == "anc_post_reset":
                    p = float(p_reset[r, apos[qb]])
                    if p > 0:
                        s = f"X_ERROR({p!r}) {qb}"
                        parts.append(s); mparts.append(s)
                    continue
                if qb in dpos:
                    g, q = float(g_data[r, dpos[qb]]), float(q_data[r, dpos[qb]])
                else:
                    g, q = float(g_anc[r, apos[qb]]), float(q_anc[r, apos[qb]])
                if g > 0 or q > 0:
                    pp = tc.layer_pauli(g, q, kw[kind])
                    if any(pp):
                        s = "PAULI_CHANNEL_1({!r},{!r},{!r}) {}".format(pp[0], pp[1], pp[2], qb)
                        parts.append(s); mparts.append(s)
            if parts:
                out += stim.Circuit("\n".join(parts))
            if mod is not None and mparts:
                mod += stim.Circuit("\n".join(mparts))
    return out, mod


def sample_trajectory_fast(circuit, shots: int, seed: int, rec_idx, anc_pos, eps0, eps1, rounds: int, chunk: int = 256):
    """All shots of one task: sampler and converter compiled once; readout flips drawn only on records whose
    eps0 or eps1 is non-zero (independent Bernoulli per record, as in sample_trajectory; a different random stream).
    Returns (det [S, n_det], obs, y [S, R, n_anc])."""
    n_anc = len(rec_idx) // rounds
    rr = np.repeat(np.arange(rounds), n_anc)
    e0, e1 = np.asarray(eps0)[rr, anc_pos], np.asarray(eps1)[rr, anc_pos]
    act = np.nonzero((e0 > 0) | (e1 > 0))[0]
    smp = circuit.compile_sampler(seed=seed)
    conv = circuit.compile_m2d_converter()
    rng = np.random.default_rng([int(seed) % (2 ** 63), 5])
    dets, obss, ys = [], [], []
    for k in range(0, shots, chunk):
        n = min(chunk, shots - k)
        meas = smp.sample(n)
        m = meas[:, rec_idx]
        if len(act):
            ma = m[:, act]
            flip = rng.random(ma.shape) < np.where(ma, e1[act], e0[act])
            m[:, act] = ma ^ flip
            meas[:, rec_idx] = m
        d_, o_ = conv.convert(measurements=meas, separate_observables=True)
        dets.append(d_); obss.append(o_); ys.append(m.reshape(n, rounds, n_anc))
    return np.concatenate(dets), np.concatenate(obss), np.concatenate(ys)
