# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "scipy"]
# ///
"""Barcode co-occurrence matrix preprocessor for the Barcode Contact Map viewer.

Writes one gzipped `.bcm` bundle per chromosome: contacts plus the hotspots and junctions documents.
Section numbers in the comments (§) refer to the design notes the tool was built from.

  uv run preprocess.py --input chr1_5kb_bin.stat [chr2_5kb_bin.stat.gz ...] --out bcm/ [--fai REF.fa.fai]
                       [--name NAME [--name NAME2 ...]] [--junction-null N]
  uv run preprocess.py --synthetic --out bcm/
  uv run preprocess.py --verify bcm/chr1.bcm [--input chr1_5kb_bin.stat]

Input: tab- or space-separated `binA_start binB_start count`, one line per bin pair, both (a, b) and (b, a)
present; plain text or gzipped.

The output is OUT/<name>.bcm per input. Name: --name (repeat once per --input, in order); else, with --fai, the
longest index name that prefixes the file name and is followed by '_', '.' or the end
(NC_000001.11_whole.stat -> NC_000001.11); else the file name up to the first '_' or '.'
(chr1_whole_5kb_bin.stat -> chr1).
"""
import argparse
import gzip
import json
import math
import re
import struct
import sys
import time
import zlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import optimize, signal, sparse, special, stats
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

MAGIC = b"BCM2"
VERSION = 2
HEADER = struct.Struct("<4s7I2d16x")  # 64 bytes, §2.2 / §10.1
U32 = struct.Struct("<I")
RESOLUTION = 50_000
FAR_BP = 1_000_000          # covFar definition, header minSepBp, near/far split (§2.1, §8.1)
TEST_MIN_SEP = 300_000      # params.minSepBp (§8.1)
LOWCOV_FACTOR = 0.1
STRIPE_FACTOR = 3
MIN_OBS = 10
MIN_OE = 3
ALPHA = 0.01
EBAL_TRIM = 0.005
PHI_TRIM = 0.001
MERGE_CHEBYSHEV = 2
POINT_SHARE = 0.9
HUB_WINDOW = 1_000_000
N_LEVELS = 7
NEGLOG10P_CAP = 1000
MAX_BIN_SIZE = 5 * 1024 * 1024
MAX_JSON_SIZE = 1024 * 1024
# §10.1 decoder sanity limits (same as the viewer): checked on the header before anything is allocated
MAX_NBINS = 1 << 20
MAX_NNZ = 1 << 29
MAX_TRAILER = 64 * 1024 * 1024  # jsonLen cap; the two documents are each <= 1 MB in practice
CALIBRATION_NOTE = ("the junction score threshold (tau = 450) and the null tail constants (params.null) were "
                    "calibrated on a single 93 Mb chromosome; for another chromosome, --junction-null N re-runs the diagnostic null "
                    "(about 90 s per 150,000 boxes) and prints the score at which expFP = 1")


def log(msg):
    print(msg, flush=True)


def sig4(x):
    """Round a float to 4 significant digits (§2.3)."""
    return float(f"{x:.4g}")


def ceil_div(a, b):
    return -(-a // b)


# ---------------------------------------------------------------- §3.1 encode

def load_tsv(path, bin_size):
    """§3.1.1: parse `a b count` lines into bin indices i, j and counts c."""
    data = Path(path).read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    try:
        arr = np.fromstring(data, dtype=np.int64, sep=" ")
    except ValueError:
        sys.exit(f"{path}: found a non-numeric line (a header?); expected `binA_start binB_start count` only")
    n_lines = data.count(b"\n") + (not data.endswith(b"\n"))
    empty = not data.strip()
    del data
    if empty:
        sys.exit(f"{path}: no data lines")
    if arr.size != 3 * n_lines:
        sys.exit(f"{path}: parsed {arr.size} numbers from {n_lines} lines; expected 3 per line")
    a, b, c = arr.reshape(-1, 3).T
    if np.any(a % bin_size) or np.any(b % bin_size):
        sys.exit(f"{path}: positions are not multiples of bin size {bin_size}")
    if np.any(c <= 0) or c.max() >= 2**31:
        sys.exit(f"{path}: counts must be in [1, 2^31)")
    return (a // bin_size).astype(np.int32), (b // bin_size).astype(np.int32), c.astype(np.int32)


def check_symmetric(i, j, c, n, bs=1):
    """§3.1.2: rows with a<b must mirror rows with a>b exactly; no duplicate cells."""
    lo, hi = i < j, i > j
    k1 = i[lo].astype(np.int64) * n + j[lo]
    k2 = j[hi].astype(np.int64) * n + i[hi]
    o1, o2 = np.argsort(k1), np.argsort(k2)
    k1, k2 = k1[o1], k2[o2]
    kd = np.sort(i[i == j])
    if np.any(k1[1:] == k1[:-1]) or np.any(kd[1:] == kd[:-1]):
        sys.exit("input has duplicate (a, b) lines")
    if k1.size != k2.size or not (np.array_equal(k1, k2) and np.array_equal(c[lo][o1], c[hi][o2])):
        only = np.setdiff1d(k1, k2)
        if only.size:
            a, b = (x * bs for x in divmod(int(only[0]), n))
            sys.exit(f"input matrix is not symmetric: line `{a} {b}` has no `{b} {a}` mirror")
        only = np.setdiff1d(k2, k1)
        if only.size:
            a, b = (x * bs for x in divmod(int(only[0]), n))
            sys.exit(f"input matrix is not symmetric: line `{b} {a}` has no `{a} {b}` mirror")
        k = int(np.flatnonzero(c[lo][o1] != c[hi][o2])[0])
        a, b = (x * bs for x in divmod(int(k1[k]), n))
        sys.exit(f"input matrix is not symmetric: lines `{a} {b}` and `{b} {a}` have different counts")


def coverages(i, j, c, n, far_bins):
    """§2.1 / §3.1.3: full-row sums covAll and raw covFar (|i-j| >= far_bins) over all lines."""
    cov_all = np.bincount(i, weights=c, minlength=n)
    far = np.abs(i - j) >= far_bins
    cov_far = np.bincount(i[far], weights=c[far], minlength=n)
    return cov_all.astype(np.int64), cov_far.astype(np.int64)


def varint_encode(values):
    """Unsigned LEB128 for a non-negative int array, vectorised (§2.2 cols)."""
    v = values.astype(np.uint64)
    nbytes = np.ones(v.size, dtype=np.int64)
    for k in range(1, 5):
        nbytes += v >= (1 << (7 * k))
    starts = np.concatenate(([0], np.cumsum(nbytes)[:-1]))
    out = np.zeros(int(nbytes.sum()), dtype=np.uint8)
    for k in range(int(nbytes.max(initial=1))):
        m = nbytes > k
        byte = (v[m] >> np.uint64(7 * k)) & np.uint64(0x7F)
        cont = (nbytes[m] > k + 1).astype(np.uint64) << np.uint64(7)
        out[starts[m] + k] = (byte | cont).astype(np.uint8)
    return out


def encode(i, j, c, n, bin_size, far_bp, cov_far):
    """§3.1.4: build the contacts payload (header + §2.2 sections 1-6) with the given (self-excluded, §8.1) covFar."""
    cov_all, _ = coverages(i, j, c, n, far_bp // bin_size)
    total = int(c.sum())
    diag = int(c[i == j].sum())
    up = i <= j
    ui, uj, uc = i[up], j[up], c[up]
    order = np.lexsort((uj, ui))
    ui, uj, uc = ui[order], uj[order], uc[order]
    nnz = ui.size

    row_ptr = np.zeros(n + 1, dtype=np.int64)
    row_ptr[1:] = np.cumsum(np.bincount(ui, minlength=n))
    delta = np.empty(nnz, dtype=np.int64)
    delta[1:] = uj[1:] - uj[:-1]
    first = np.zeros(nnz, dtype=bool)
    first[row_ptr[:-1][np.diff(row_ptr) > 0]] = True
    delta[first] = uj[first] - ui[first]
    cols = varint_encode(delta)

    counts = np.minimum(uc, 255).astype(np.uint8)
    overflow = uc[uc >= 255]
    pad = (-nnz) % 4

    for name, arr in (("rowPtr", row_ptr), ("covAll", cov_all), ("covFar", cov_far), ("overflow", overflow)):
        if arr.size and (arr.max() >= 2**32 or arr.min() < 0):
            sys.exit(f"{name} out of u32 range")

    header = HEADER.pack(MAGIC, VERSION, bin_size, n, nnz, overflow.size, cols.size, far_bp,
                         float(total), float(diag))
    payload = b"".join([
        header,
        row_ptr.astype("<u4").tobytes(),
        cov_all.astype("<u4").tobytes(),
        cov_far.astype("<u4").tobytes(),
        counts.tobytes(), b"\0" * pad,
        overflow.astype("<u4").tobytes(),
        cols.tobytes(),
    ])
    info = dict(nnz=nnz, totalCount=total, diagCount=diag, nOverflow=int(overflow.size), colBytes=int(cols.size))
    return payload, info


def bundle_json(name, chrom_length, hot_doc, jun_doc):
    """§10.1 trailer document (compact UTF-8 JSON; NaN/Infinity are rejected so JSON.parse can read it)."""
    doc = dict(format="bcm", version=VERSION, chrom=name, chromLength=chrom_length, createdBy="preprocess.py",
               hotspots=hot_doc, junctions=jun_doc)
    return json.dumps(doc, separators=(",", ":"), allow_nan=False).encode("utf-8")


def bundle_payload(contacts, col_bytes, meta):
    """§10.1: contacts payload (header + sections 1-6), zero pad cols to 4 bytes, u32 jsonLen, JSON bytes."""
    return b"".join([contacts, b"\0" * ((-col_bytes) % 4), U32.pack(len(meta)), meta])


# ------------------------------------------------- §2.1 / §8.1 background model

def block_sum(v, fac, n_out):
    """Sum consecutive groups of `fac` bins (coarse coverage, §2.1)."""
    return np.bincount(np.arange(v.size) // fac, weights=v, minlength=n_out)


def coverage_model(cov_far, partial_last):
    """§2.1 + §8.1: f = covFar / median(covFar > 0); lowcov f < 0.1, stripe f >= 3.

    A partially filled last coarse bin is excluded from valid (its f is not comparable);
    it is not flagged lowcov or stripe.
    """
    pos = cov_far[cov_far > 0]
    f = cov_far / (float(np.median(pos)) if pos.size else 1.0)
    lowcov = f < LOWCOV_FACTOR
    stripe = f >= STRIPE_FACTOR
    valid = ~lowcov & ~stripe
    if partial_last:
        valid[-1] = False
    return f, lowcov, stripe, valid


def trimmed_mean(x, trim=EBAL_TRIM):
    """§8.1: drop the largest ceil(trim * count) values, average the rest (0 if nothing is left)."""
    keep = x.size - math.ceil(trim * x.size)
    if keep <= 0:
        return 0.0
    return float(np.partition(x, keep - 1)[:keep].mean())


def expected_by_distance(C, f, valid, far_d):
    """§8.1: Ebal(d) per-diagonal trimmed mean for d < far_d, one pooled trimmed far mean for d >= far_d."""
    n = C.shape[0]
    fv = np.where(valid, f, 1.0)
    ebal = np.zeros(n)
    for d in range(min(far_d, n)):
        a = np.flatnonzero(valid[:n - d] & valid[d:])
        ebal[d] = trimmed_mean(C[a, a + d] / (fv[a] * fv[a + d]))
    ti, tj = np.triu_indices(n, k=far_d)
    ok = valid[ti] & valid[tj]
    ti, tj = ti[ok], tj[ok]
    far = trimmed_mean(C[ti, tj] / (f[ti] * f[tj]))
    ebal[far_d:] = far
    return ebal, far


def drop_top_oe(O, E, n_drop):
    """Boolean keep-mask removing the n_drop cells with the largest O/E."""
    keep = np.ones(O.size, dtype=bool)
    n_drop = min(n_drop, O.size)
    if n_drop > 0:
        keep[np.argpartition(-(O / E), n_drop - 1)[:n_drop]] = False
    return keep


def nb_mle_phi(O, E):
    """§8.1: NB maximum-likelihood phi = 1/k with offsets E, after dropping the top 0.1% by O/E."""
    m = E > 0
    O, E = O[m], E[m]
    if O.size == 0:
        return 0.0
    keep = drop_top_oe(O, E, math.ceil(PHI_TRIM * O.size))
    O, E = O[keep], E[keep]
    logE = np.log(E)

    def nll(logk):
        k = math.exp(logk)
        return -np.sum(special.gammaln(O + k) - special.gammaln(k) - k * np.log1p(E / k)
                       + O * (logE - np.log(k + E)))
    return math.exp(-optimize.minimize_scalar(nll, bounds=(-6, 12), method="bounded").x)


def neg_log10_p(O, E, phi):
    """§3.2.5: -log10 P(X >= O) under NB(mean E, size 1/phi) (phi may be per cell)."""
    k = 1.0 / np.asarray(phi)
    return -stats.nbinom.logsf(O - 1, k, k / (k + E)) / math.log(10)


def phi_by_level(ui, uj, uc, cov_far5, n, bin_size, ebal_far50k):
    """§8.1 params.phiByLevel: trimmed ratio-of-sums sum((O-E)^2-O)/sum(E^2) over far valid cells, per level.

    Uses the viewer's rules: f from covFar at that level, Ebal_far = ebalFar50k * (bs/50000)^2.
    Zero cells enter through closed-form sums of E^2 over all valid far pairs.
    """
    out = []
    for L in range(N_LEVELS):
        fac = 1 << L
        bs, nL = bin_size * fac, ceil_div(n, fac)
        f, _, _, valid = coverage_model(block_sum(cov_far5, fac, nL), n % fac != 0)
        efar = ebal_far50k * (bs / RESOLUTION) ** 2
        far_d = ceil_div(FAR_BP, bs)
        I, J = ui >> L, uj >> L
        m = (J - I >= far_d) & valid[I] & valid[J]
        key, inv = np.unique(I[m].astype(np.int64) * nL + J[m], return_inverse=True)
        O = np.bincount(inv, weights=uc[m])
        E = efar * f[key // nL] * f[key % nL]

        g = np.where(valid, f * f, 0.0)
        suf_g = np.concatenate((np.cumsum(g[::-1])[::-1], [0.0]))
        suf_v = np.concatenate((np.cumsum(valid[::-1])[::-1], [0]))
        nxt = np.minimum(np.arange(nL) + far_d, nL)
        sum_e2 = efar ** 2 * float(np.sum(g * suf_g[nxt]))
        n_pairs = int(np.sum(valid * suf_v[nxt]))

        keep = drop_top_oe(O, E, math.ceil(PHI_TRIM * O.size))  # trim nonzero cells only
        den = sum_e2 - float(np.sum(E[~keep] ** 2))
        num = float(np.sum((O * O - 2 * O * E - O)[keep])) + den
        out.append(max(0.0, num / den) if den > 0 else 0.0)
    return out


# ------------------------------------------------------- §3.2 hotspots at 50 kb

def test_cells(C, cov_far, partial_last, min_d, far_d):
    """§3.2 steps 2-6 with §8.1 changes: fit f, Ebal, phiNear/phiFar and NB-test every eligible cell."""
    n = C.shape[0]
    f, lowcov, stripe, valid = coverage_model(cov_far, partial_last)
    ebal, ebal_far = expected_by_distance(C, f, valid, far_d)
    ti, tj = np.triu_indices(n, k=min_d)
    ok = valid[ti] & valid[tj]
    ti, tj = ti[ok], tj[ok]
    O = C[ti, tj]
    E = ebal[tj - ti] * f[ti] * f[tj]
    near = (tj - ti) < far_d
    phi_near, phi_far = nb_mle_phi(O[near], E[near]), nb_mle_phi(O[~near], E[~near])
    phi = np.where(near, phi_near, phi_far)
    n_tests = int(ti.size)
    with np.errstate(divide="ignore", invalid="ignore"):
        oe = np.where(E > 0, O / E, np.inf)
    cand = (O >= MIN_OBS) & (oe >= MIN_OE)
    nlp = np.zeros(n_tests)
    nlp[cand] = neg_log10_p(O[cand], E[cand], phi[cand])
    sig = cand & (nlp >= -math.log10(ALPHA / max(n_tests, 1)))
    return SimpleNamespace(f=f, lowcov=lowcov, stripe=stripe, valid=valid, ebal=ebal, ebal_far=ebal_far,
                           ti=ti, tj=tj, O=O, E=E, nlp=nlp, cand=cand, sig=sig, near=near,
                           phi_near=phi_near, phi_far=phi_far, n_tests=n_tests)


def exclude_significant(ui, uj, uc, cov_far5, si, sj, fac, far_bins):
    """§8.1 self-exclusion: remove far 5 kb cells inside significant 50 kb cells from covFar (row and column)."""
    nr = ceil_div(cov_far5.size, fac)
    hit = np.zeros((nr, nr), dtype=bool)
    hit[si, sj] = True
    m = hit[ui // fac, uj // fac] & (uj - ui >= far_bins)
    n = cov_far5.size
    sub = np.bincount(ui[m], weights=uc[m], minlength=n) + np.bincount(uj[m], weights=uc[m], minlength=n)
    return cov_far5 - sub.astype(np.int64)


def merge_cells(si, sj):
    """§3.2.7: connected components of significant cells with Chebyshev distance <= 2."""
    if si.size == 0:
        return np.zeros(0, dtype=np.int64)
    pairs = cKDTree(np.column_stack((si, sj))).query_pairs(MERGE_CHEBYSHEV, p=np.inf, output_type="ndarray")
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(si.size, si.size))
    return connected_components(g, directed=False)[1]


def build_regions(t, data_len):
    """§3.2.8: one record per merged component, sorted by excess, ranked from 1."""
    s = t.sig
    si, sj, sO, sE, snlp = t.ti[s], t.tj[s], t.O[s], t.E[s], t.nlp[s]
    labels = merge_cells(si, sj)
    out = []
    for lab in range(labels.max(initial=-1) + 1):
        m = labels == lab
        ci, cj, cO, cE, cp = si[m], sj[m], sO[m], sE[m], snlp[m]
        pk = np.lexsort((-cO, -cp))[0]
        obs, exp = float(cO.sum()), float(cE.sum())
        out.append(dict(
            x0=int(cj.min() * RESOLUTION), x1=int(min((cj.max() + 1) * RESOLUTION, data_len)),
            y0=int(ci.min() * RESOLUTION), y1=int(min((ci.max() + 1) * RESOLUTION, data_len)),
            peakX=int(cj[pk] * RESOLUTION), peakY=int(ci[pk] * RESOLUTION),
            obs=int(round(obs)), exp=sig4(exp), excess=obs - exp,
            peakOE=sig4(cO[pk] / cE[pk]) if cE[pk] > 0 else None,
            peakNegLog10P=sig4(min(cp[pk], NEGLOG10P_CAP)), nCells=int(m.sum()),
        ))
    out.sort(key=lambda h: -h["excess"])
    return [{"rank": r, **h, "excess": sig4(h["excess"])} for r, h in enumerate(out, 1)]


def shape_fields(hs, i, j, c, bin_size):
    """§8.1: nz5k, maxRowShare5k, maxColShare5k and shape from the 5 kb cells inside each bbox.

    i, j, c are full (mirrored) 5 kb lines.
    """
    order = np.argsort(i, kind="stable")
    fi, fj, fc = i[order], j[order], c[order]
    ptr = np.searchsorted(fi, np.arange(int(fi.max(initial=0)) + 2))
    for h in hs:
        r0, r1 = h["y0"] // bin_size, ceil_div(h["y1"], bin_size)
        c0, c1 = h["x0"] // bin_size, ceil_div(h["x1"], bin_size)
        a, b = ptr[min(r0, ptr.size - 1)], ptr[min(r1, ptr.size - 1)]
        m = (fj[a:b] >= c0) & (fj[a:b] < c1)
        rr, cc, vv = fi[a:b][m] - r0, fj[a:b][m] - c0, fc[a:b][m]
        tot = float(vv.sum())
        row_share = np.bincount(rr, weights=vv).max(initial=0) / tot if tot else 0.0
        col_share = np.bincount(cc, weights=vv).max(initial=0) / tot if tot else 0.0
        h["nz5k"] = int(vv.size)
        h["maxRowShare5k"] = sig4(row_share)
        h["maxColShare5k"] = sig4(col_share)
        h["shape"] = "point" if max(row_share, col_share) > POINT_SHARE else "extended"


def hub_degrees(hs):
    """§8.1 hubDegree: other hotspots for which some 1 Mb window meets an interval of both (either axis)."""
    if not hs:
        return np.zeros(0, dtype=np.int64)
    iv = np.array([[[h["y0"], h["y1"]], [h["x0"], h["x1"]]] for h in hs], dtype=np.int64)
    a, b = iv[:, None, :, None, :], iv[None, :, None, :, :]
    gap = np.maximum(a[..., 0], b[..., 0]) - np.minimum(a[..., 1], b[..., 1])
    share = (gap < HUB_WINDOW).any(axis=(2, 3))
    np.fill_diagonal(share, False)
    return share.sum(axis=1)


def runs(mask):
    """Start/end (exclusive) indices of runs of True in a 1-D boolean array."""
    m = np.concatenate(([False], mask, [False])).astype(np.int8)
    d = np.diff(m)
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def detect_hotspots(i, j, c, n, bin_size, far_bp, test_min_sep, chrom_length):
    """§3.2 + §8.1: two-pass NB calling at 50 kb. Returns JSON parts and the self-excluded 5 kb covFar.

    Coordinates are clipped to chrom_length (= n * bin_size unless --fai / --chrom-length gave a shorter last bin).
    """
    fac = RESOLUTION // bin_size
    nr, partial = ceil_div(n, fac), n % fac != 0
    data_len = min(n * bin_size, chrom_length)
    C = np.bincount((i // fac).astype(np.int64) * nr + j // fac, weights=c, minlength=nr * nr).reshape(nr, nr)
    up = i <= j
    ui, uj, uc = i[up], j[up], c[up]
    far_bins = far_bp // bin_size
    _, cov_far5 = coverages(i, j, c, n, far_bins)
    min_d, far_d = ceil_div(test_min_sep, RESOLUTION), ceil_div(far_bp, RESOLUTION)

    t = None
    for p in (1, 2):
        if t is not None:
            cov_far5 = exclude_significant(ui, uj, uc, cov_far5, t.ti[t.sig], t.tj[t.sig], fac, far_bins)
        t = test_cells(C, block_sum(cov_far5, fac, nr), partial, min_d, far_d)
        log(f"pass {p}: res {RESOLUTION} n {nr} nTests {t.n_tests} (near {int(t.near.sum())}) "
            f"phiNear {t.phi_near:.4f} phiFar {t.phi_far:.4f} ebalFar {t.ebal_far:.4f} "
            f"candidates {int(t.cand.sum())} significant {int(t.sig.sum())} (near {int((t.sig & t.near).sum())})")

    hotspots = build_regions(t, data_len)
    shape_fields(hotspots, i, j, c, bin_size)
    for h, k in zip(hotspots, hub_degrees(hotspots)):
        h["hubDegree"] = int(k)

    s0, s1 = runs(t.stripe)
    stripes = [dict(start=int(a * RESOLUTION), end=int(min(b * RESOLUTION, data_len)),
                    maxFactor=sig4(t.f[a:b].max())) for a, b in zip(s0, s1)]
    l0, l1 = runs(t.lowcov)
    lowcov = [[int(a * RESOLUTION), int(min(b * RESOLUTION, data_len))] for a, b in zip(l0, l1)]
    pbl = phi_by_level(ui, uj, uc, cov_far5, n, bin_size, sig4(t.ebal_far))
    log(f"hotspots: {len(hotspots)} regions, {sum(h['excess'] >= 100 for h in hotspots)} with excess >= 100, "
        f"{sum(h['excess'] >= 500 for h in hotspots)} >= 500, {sum(h['shape'] == 'point' for h in hotspots)} point; "
        f"{len(stripes)} stripes; {len(lowcov)} lowcov runs ({int(t.lowcov.sum())}/{nr} bins)")
    log("phiByLevel (5..320 kb): " + " ".join(f"{x:.4g}" for x in pbl))
    params = dict(resolution=RESOLUTION, minSepBp=test_min_sep, lowcovFactor=LOWCOV_FACTOR,
                  stripeFactor=STRIPE_FACTOR, minObs=MIN_OBS, minOE=MIN_OE, alpha=ALPHA, nTests=t.n_tests,
                  phi=sig4(t.phi_far), phiNear=sig4(t.phi_near), phiFar=sig4(t.phi_far),
                  phiByLevel=[sig4(x) for x in pbl], ebalFar50k=sig4(t.ebal_far),
                  mergeChebyshev=MERGE_CHEBYSHEV)
    return params, hotspots, stripes, lowcov, cov_far5


# ----------------------------------------------------------- §3.3 synthetic

SYN_N, SYN_BIN, SYN_NAME = 2000, 5000, "synthetic"
SYN_EVENTS = [(600, 900), (1200, 1700)]  # (row, col) bins: inversion corner, translocation
SYN_STRIPE = (1500, 1504)
SYN_GAP = (100, 140)
# DESIGN §3.3 says 300 and 200. At that strength the event rows/cols reach f >= 3 against the
# 0.02 background, get classified as stripes and are never tested, so the fixture had no hotspots.
SYN_INV_AMP, SYN_TRANS_AMP = 30, 20


def synthetic_tsv():
    """§3.3: seed-1 Poisson fixture, mirrored to full symmetric (i, j, c) lines."""
    rng = np.random.default_rng(1)
    n = SYN_N
    r, q = np.indices((n, n))
    d = q - r
    lam = np.full((n, n), 0.02)
    near = (d >= 0) & (d <= 30)
    lam[near] += 1000 * np.exp(-d[near] / 2)
    u, v = np.indices((21, 21))
    b1, b2 = SYN_EVENTS[0]
    np.add.at(lam, (b1 - u, b2 - v), SYN_INV_AMP * np.exp(-(u + v) / 3))
    np.add.at(lam, (b1 + u, b2 + v), SYN_INV_AMP * np.exp(-(u + v) / 3))
    t1, t2 = SYN_EVENTS[1]
    np.add.at(lam, (t1 + u, t2 - v), SYN_TRANS_AMP * np.exp(-(u + v) / 3))
    s = np.zeros(n, dtype=bool)
    s[SYN_STRIPE[0]:SYN_STRIPE[1]] = True
    sm = s[:, None] | s[None, :]
    lam[sm] = lam[sm] * 8 + 0.2
    counts = rng.poisson(lam)
    g = np.zeros(n, dtype=bool)
    g[SYN_GAP[0]:SYN_GAP[1]] = True
    counts[g[:, None] | g[None, :]] = 0
    counts = np.triu(counts)
    ui, uj = np.nonzero(counts)
    uc = counts[ui, uj]
    off = ui != uj
    i = np.concatenate((ui, uj[off]))
    j = np.concatenate((uj, ui[off]))
    c = np.concatenate((uc, uc[off]))
    return i.astype(np.int32), j.astype(np.int32), c.astype(np.int32), n


# ---------------------------------------------------------- §9 junctions at 5 kb

JUN_L = 20                  # template length (bins, 100 kb)
JUN_G = 4                   # seed box side (bins)
JUN_R = 2                   # proxy margin around a seed box (bins)
SEED_OBS, SEED_OE = 10, 5
TAU, TAU_SCREEN, TAU2 = 450.0, 100.0, 50.0
CORE_MIN, CORE_FRAC = 25.0, 0.1
INNER_D, INNER_SUP = 4, 0.15
LINE_POINT = 0.33
LINE_GAP = 4                # line profiles skip cells within 4 bins of the vertex on the other axis
NMS_HALF = 20               # half-bins (50 kb)
POLISH_HALF = 4
UNC_DROP = 3.0
CLIMB_STEPS = 6
FIT_ITERS = 8
PHI_SAMPLE = 400_000
PHI_FLOOR = 1e-6
JUN_NULL = dict(u=100.0, beta=65.0, pU=3.45e-4, basin=0.096)
FP_FLOOR = 1e-12
QUADS = ("--", "++", "-+", "+-")
PAIRS = ((0, 1), (2, 3))    # {--, ++}, {-+, +-}
SINGLE_EVENT = {"--": "inversion-left", "++": "inversion-right", "-+": "deletion", "+-": "tandem-dup"}
JUN_TYPES = ("point", "junction", "butterfly-inversion", "butterfly-deldup", "block")


class JunctionModel:
    """§9.2: 5 kb background (f, valid, Ebal5, phiNear/phiFar), template and per-quadrant NB fits.

    Vertices are half-bin pairs (a2, b2); all methods are vectorised over arrays of vertices.
    """

    def __init__(self, ui, uj, uc, cov_far5, n, bs, hot_params):
        self.n, self.bs = n, bs
        self.dmin = ceil_div(hot_params["minSepBp"], bs)
        self.far_d = ceil_div(FAR_BP, bs)
        U = csr_matrix((uc.astype(np.float64), (ui, uj)), shape=(n, n))
        U.sum_duplicates()
        self.U = U
        self.M = (U + sparse.triu(U, 1).T).tocsr()
        self.f, _, _, self.valid = coverage_model(cov_far5.astype(np.float64), partial_last=False)
        self.ebal_far = hot_params["ebalFar50k"] * (bs / RESOLUTION) ** 2
        self.ebal = self.expected_5k(U.tocoo())
        self.t = np.zeros(JUN_L)  # chromosomes shorter than L + 1 bins: missing distances get t = 0
        e = self.ebal[1:JUN_L + 1]
        self.t[:e.size] = np.maximum(e - self.ebal_far, 0)
        self.phi_near, self.phi_far = self.phi_bands()
        self.k_near, self.k_far = 1 / self.phi_near, 1 / self.phi_far
        u, v = np.nonzero(np.add.outer(np.arange(JUN_L), np.arange(JUN_L)) <= JUN_L - 1)
        self.u, self.v, self.tt = u, v, self.t[u + v]
        self.inner = (u + v + 1) <= INNER_D
        self.core = ((u >= 1) & (v >= 1))[None]

    def expected_5k(self, U):
        """§9.2 Ebal5: per-diagonal 99.5% trimmed mean (zeros included) from sparse cells; ebalFar5 beyond 1 Mb."""
        r, c, v = U.row, U.col, U.data
        d = c - r
        f, valid, n = self.f, self.valid, self.n
        ok = valid[r] & valid[c] & (d < self.far_d)
        bal, d = v[ok] / (f[r[ok]] * f[c[ok]]), d[ok]
        order = np.argsort(d, kind="stable")
        bal, bounds = bal[order], np.searchsorted(d[order], np.arange(self.far_d + 1))
        ebal = np.full(n, self.ebal_far)
        for k in range(min(self.far_d, n)):
            cnt = int(np.sum(valid[:n - k] & valid[k:]))
            drop = math.ceil(EBAL_TRIM * cnt)
            if cnt - drop <= 0:
                ebal[k] = 0
                continue
            x = bal[bounds[k]:bounds[k + 1]]
            if drop:
                x = np.sort(x)[:-drop] if drop < x.size else np.zeros(0)
            ebal[k] = x.sum() / (cnt - drop)
        return ebal

    def lookup(self, r, c):
        if r.size == 0:  # scipy returns a 1 x 0 matrix for empty fancy indices
            return np.zeros(r.shape)
        return np.asarray(self.M[r.ravel(), c.ravel()]).reshape(r.shape)

    def phi_bands(self):
        """§9.2 dispersion: NB MLE on 400,000 random valid cells per band (near first), one shared rng."""
        rng, n, out = np.random.default_rng(0), self.n, []
        for lo, hi in ((self.dmin, self.far_d), (self.far_d, n)):
            if hi <= lo or n == 0:  # chromosome too short for this band (the far band is drawn last)
                out.append(PHI_FLOOR)
                continue
            r = rng.integers(0, n, 8 * PHI_SAMPLE)
            c = r + rng.integers(lo, hi, 8 * PHI_SAMPLE)
            r, c = r[c < n], c[c < n]
            ok = self.valid[r] & self.valid[c]
            r, c = r[ok][:PHI_SAMPLE], c[ok][:PHI_SAMPLE]
            O, E = self.lookup(r, c), self.ebal[c - r] * self.f[r] * self.f[c]
            out.append(max(nb_mle_phi(O, E), PHI_FLOOR) if O.size else PHI_FLOOR)
        return out

    def window(self, r0, r1, c0, c1):
        """Dense O, Exp, validity and |c - r| over rows [r0, r1) x cols [c0, c1); out-of-range cells invalid."""
        n = self.n
        rr, cc = np.arange(r0, r1), np.arange(c0, c1)
        inr, inc = (rr >= 0) & (rr < n), (cc >= 0) & (cc < n)
        O = np.zeros((rr.size, cc.size))
        a0, a1, b0, b1 = max(r0, 0), min(r1, n), max(c0, 0), min(c1, n)
        if a1 > a0 and b1 > b0:
            O[np.ix_(inr, inc)] = self.M[a0:a1, b0:b1].toarray()
        rc, ccl = np.clip(rr, 0, n - 1), np.clip(cc, 0, n - 1)
        d = np.abs(ccl[None, :] - rc[:, None])
        E = self.ebal[d] * self.f[rc][:, None] * self.f[ccl][None, :]
        V = (inr & self.valid[rc])[:, None] & (inc & self.valid[ccl])[None, :]
        return O, E, V, d

    def cells(self, a2, b2, q):
        """§9.2 quadrant cells for vertices (a2, b2): O, Exp, template T, usable mask W, size k, tcov."""
        n, odd_a, odd_b = self.n, a2 % 2 == 1, b2 % 2 == 1
        A = a2 // 2 + (odd_a & (q[0] == "-"))
        B = b2 // 2 + (odd_b & (q[1] == "-"))
        r = A[:, None] + (self.u if q[0] == "+" else -1 - self.u)[None]
        c = B[:, None] + (self.v if q[1] == "+" else -1 - self.v)[None]
        # vertices with B - A < dMinBp are never scored, whatever the quadrant
        ok = (r >= 0) & (c < n) & (c - r >= self.dmin) & (b2 - a2 >= 2 * self.dmin)[:, None]
        if q[0] == "-":  # both odd: the shared centre cell belongs to ++ and +- only
            ok &= ~((odd_a & odd_b)[:, None] & ((self.u == 0) & (self.v == 0))[None])
        r, c = np.clip(r, 0, n - 1), np.clip(c, 0, n - 1)
        ff = self.f[r] * self.f[c]
        E = self.ebal[np.clip(c - r, 0, n - 1)] * ff
        W = ok & self.valid[r] & self.valid[c] & (E > 0)  # Exp = 0 cells (sparse inputs) are unusable
        E = np.where(W, E, 1.0)
        T = np.where(W, self.tt[None] * ff, 0.0)
        O = np.where(W, self.lookup(r, c), 0.0)
        K = np.where(c - r >= self.far_d, self.k_far, self.k_near)
        tcov = np.sum(np.where(W, self.tt[None], 0), axis=1) / self.tt.sum()
        return O, E, T, W, K, tcov

    @staticmethod
    def fit(O, E, T, W, K):
        """§9.2: amplitude s >= 0 of O ~ NB(Exp + s T, k) by Fisher scoring; returns s and llr (floored at 0)."""
        W = W & (E > 0)
        E = np.where(W, E, 1.0)
        num = np.sum(np.where(W, T * (O - E) / E, 0), axis=1)
        den = np.sum(np.where(W, T * T / E, 0), axis=1)
        s = np.where(den > 0, np.maximum(num, 0) / np.maximum(den, 1e-300), 0.0)
        for _ in range(FIT_ITERS):
            mu = E + s[:, None] * T
            w = np.where(W, K / (mu * (K + mu)), 0)
            g, info = np.sum(T * (O - mu) * w, axis=1), np.sum(T * T * w, axis=1)
            sn = s + np.where(info > 0, g / np.maximum(info, 1e-300), 0)
            s = np.maximum(np.where(s > 0, np.clip(sn, s / 8, s * 8), sn), 0)
        mu = E + s[:, None] * T
        llr = np.sum(np.where(W, K * np.log((K + E) / (K + mu)) + O * np.log(mu * (K + E) / (E * (K + mu))), 0),
                     axis=1)
        return s, np.maximum(llr, 0)

    def quad_stats(self, a2, b2, chunk=4000):
        """§9.2 per-quadrant s, llr, core, obs, exp, inner; arrays [N, 4] in QUADS order."""
        N = a2.size
        out = {k: np.zeros((N, 4)) for k in ("s", "llr", "core", "obs", "exp", "inner")}
        for k0 in range(0, N, chunk):
            sl = slice(k0, k0 + chunk)
            for jq, q in enumerate(QUADS):
                O, E, T, W, K, tcov = self.cells(a2[sl], b2[sl], q)
                s, llr = self.fit(O, E, T, W, K)
                _, core = self.fit(O, E, T, W & self.core, K)
                ok = tcov >= 0.5
                X = np.where(W, np.maximum(O - E, 0), 0)
                tot = X.sum(axis=1)
                out["s"][sl, jq] = np.where(ok, s, 0)
                out["llr"][sl, jq] = np.where(ok, llr, 0)
                out["core"][sl, jq] = np.where(ok, core, 0)
                out["obs"][sl, jq] = O.sum(axis=1)
                out["exp"][sl, jq] = np.sum(np.where(W, E, 0), axis=1)
                out["inner"][sl, jq] = np.where(tot > 0, np.sum(np.where(self.inner[None], X, 0), axis=1)
                                                / np.maximum(tot, 1e-12), 0)
        return out

    @staticmethod
    def supported(st):
        """§9.2 support rule: llr >= TAU2, core >= max(25, 0.1 llr), inner >= 0.15."""
        llr = st["llr"]
        return (llr >= TAU2) & (st["core"] >= np.maximum(CORE_MIN, CORE_FRAC * llr)) & (st["inner"] >= INNER_SUP)

    @classmethod
    def choose(cls, st):
        """§9.2 corner score and detected-quadrant mask: best single, {--,++} or {-+,+-}; else best raw llr."""
        sup = cls.supported(st)
        le = np.where(sup, st["llr"], 0.0)
        N = le.shape[0]
        j1 = le.argmax(axis=1)
        score = le[np.arange(N), j1]
        mask = np.zeros((N, 4), bool)
        mask[np.arange(N), j1] = score > 0
        for a, b in PAIRS:
            pair = le[:, a] + le[:, b]
            better = sup[:, a] & sup[:, b] & (pair > score)
            score = np.where(better, pair, score)
            mask[better] = False
            mask[better, a] = True
            mask[better, b] = True
        none = ~mask.any(axis=1)
        return np.where(none, st["llr"].max(axis=1), score), mask

    def score(self, a2, b2):
        return self.choose(self.quad_stats(a2, b2))[0]

    def seed_boxes(self):
        """§9.3.1: 4 x 4 boxes over usable upper cells with O >= 10 and O >= 5 Ebox."""
        U = self.U.tocoo()
        r, q, v = U.row, U.col, U.data
        ok = (q - r >= self.dmin) & self.valid[r] & self.valid[q]
        r, q, v = r[ok], q[ok], v[ok]
        nb = self.n // JUN_G + 1
        uk, inv = np.unique((r // JUN_G).astype(np.int64) * nb + q // JUN_G, return_inverse=True)
        O = np.bincount(inv, weights=v)
        I, J = uk // nb, uk % nb
        F = np.bincount(np.arange(self.n) // JUN_G, weights=self.f * self.valid, minlength=nb)
        Ebox = self.ebal[np.clip((J - I) * JUN_G, 0, self.n - 1)] * F[I] * F[J]
        s = (O >= SEED_OBS) & (O >= SEED_OE * Ebox)
        return I[s].astype(np.int64), J[s].astype(np.int64)

    def eligible_boxes(self):
        """§9.4 nBoxes: 4 x 4 boxes with both box starts valid and box-diagonal offset >= DMIN."""
        nb = self.n // JUN_G
        vb = self.valid[np.minimum(np.arange(nb) * JUN_G, self.n - 1)]
        cnt = np.concatenate((np.cumsum(vb[::-1])[::-1], [0]))
        nxt = np.minimum(np.arange(nb) + ceil_div(self.dmin, JUN_G), nb)
        return int(cnt[nxt][vb].sum())

    def proxy_vertices(self, I, J, batch=1024):
        """§9.3.2: start vertex per box = argmax over a 17 x 17 half-bin grid of the linear proxy Z."""
        L, G, R = JUN_L, JUN_G, JUN_R
        win = G + 2 * R + 2 * L + 2
        r0, c0 = I * G - R - L - 1, J * G - R - L - 1
        kern = {}
        for q in QUADS:
            du = self.u if q[0] == "+" else -1 - self.u
            dv = self.v if q[1] == "+" else -1 - self.v
            kern[q] = np.zeros((2 * L + 1, 2 * L + 1))
            kern[q][du + L, dv + L] = self.tt
        best = np.zeros((I.size, 2), np.int64)
        for s0 in range(0, I.size, batch):
            sl = slice(s0, s0 + batch)
            Y = np.zeros((I[sl].size, win, win))
            for k, (a, b) in enumerate(zip(r0[sl], c0[sl])):
                O, E, V, d = self.window(a, a + win, b, b + win)
                Y[k] = np.where(V & (d >= self.dmin), O - E, 0.0)
            X = {q: np.maximum(signal.fftconvolve(Y, kern[q][None, ::-1, ::-1], mode="valid", axes=(1, 2)), 0)
                 for q in QUADS}
            H = 2 * (X["++"].shape[1] - 1) - 1
            hh = np.arange(H)
            Z = {}
            for q in QUADS:
                ri = hh // 2 + ((hh % 2 == 1) & (q[0] == "-"))
                ci = hh // 2 + ((hh % 2 == 1) & (q[1] == "-"))
                Z[q] = X[q][:, ri][:, :, ci]
            tot = np.maximum(np.maximum(np.maximum(Z["--"], Z["++"]), np.maximum(Z["-+"], Z["+-"])),
                             np.maximum(Z["--"] + Z["++"], Z["-+"] + Z["+-"]))
            k = tot.reshape(tot.shape[0], -1).argmax(axis=1)
            best[sl, 0] = 2 * (r0[sl] + L) + k // H
            best[sl, 1] = 2 * (c0[sl] + L) + k % H
        return best

    def climb(self, a2, b2, cur):
        """§9.3.4: up to 6 moves to the best of the 8 half-bin neighbours while it improves by > 1e-6."""
        act = np.ones(a2.size, bool)
        da, db = (x.ravel() for x in np.meshgrid([-1, 0, 1], [-1, 0, 1], indexing="ij"))
        nz = (da != 0) | (db != 0)
        da, db = da[nz], db[nz]
        for _ in range(CLIMB_STEPS):
            idx = np.flatnonzero(act)
            if idx.size == 0:
                break
            sc = self.score((a2[idx, None] + da[None]).ravel(), (b2[idx, None] + db[None]).ravel())
            sc = sc.reshape(idx.size, 8)
            k = sc.argmax(axis=1)
            sv = sc[np.arange(idx.size), k]
            up = sv > cur[idx] + 1e-6
            a2[idx[up]] += da[k[up]]
            b2[idx[up]] += db[k[up]]
            cur[idx[up]] = sv[up]
            act[idx[~up]] = False
        return a2, b2, cur

    def polish(self, a2, b2):
        """§9.3.6: best point within +-4 half-bins; vertexUnc = max Chebyshev distance to points within 3.0.

        Exact ties (relative 1e-9) go to the point with fewer odd half-bin coordinates: on a '-' axis
        a2 = 2A - 1 and 2A give the same cells, and the vertex belongs on the bin boundary 2A.
        """
        if a2.size == 0:
            return a2.copy(), b2.copy(), np.zeros(0), np.zeros(0, np.int64)
        d = np.arange(-POLISH_HALF, POLISH_HALF + 1)
        da, db = (x.ravel() for x in np.meshgrid(d, d, indexing="ij"))
        A2, B2 = a2[:, None] + da[None], b2[:, None] + db[None]
        sc = self.score(A2.ravel(), B2.ravel()).reshape(A2.shape)
        top = sc.max(axis=1, keepdims=True)
        tied = sc >= top - 1e-9 * np.abs(top)
        k = np.where(tied, A2 % 2 + B2 % 2, 3).argmin(axis=1)
        best = sc[np.arange(a2.size), k]
        near = sc >= best[:, None] - UNC_DROP
        dist = np.maximum(np.abs(da - da[k][:, None]), np.abs(db - db[k][:, None]))
        unc = np.max(np.where(near, dist, 0), axis=1)
        return a2 + da[k], b2 + db[k], best, unc

    def line_share(self, a2, b2):
        """§9.3.8: (rowShare, colShare) per vertex from X = max(O - Exp, 0) in a +-L window."""
        L, out = JUN_L, np.zeros((a2.size, 2))
        off = np.abs(np.arange(-L, L + 1)) >= LINE_GAP
        for k in range(a2.size):
            a, b = a2[k] // 2, b2[k] // 2
            O, E, V, d = self.window(a - L, a + L + 1, b - L, b + L + 1)
            X = np.where(V & (d >= self.dmin), np.maximum(O - E, 0), 0)
            tot = X.sum()
            if tot <= 0:
                continue
            rows = np.where(off[None, :], X, 0).sum(axis=1)
            cols = np.where(off[:, None], X, 0).sum(axis=0)
            out[k] = ((rows[:-1] + rows[1:]).max() / tot, (cols[:-1] + cols[1:]).max() / tot)
        return out


def nms(a2, b2, score):
    """§9.3.5: greedy by score (stable); drop vertices within Chebyshev NMS_HALF of a kept one."""
    keep = []
    for k in np.argsort(-score, kind="stable"):
        if all(max(abs(a2[k] - a2[j]), abs(b2[k] - b2[j])) > NMS_HALF for j in keep):
            keep.append(k)
    return np.array(keep, dtype=np.int64)


def classify(mask, share):
    """§9.3.9: point (lineShare >= 0.33), block (no supported quadrant), butterfly-* or junction."""
    out = []
    for m, sh in zip(mask, share):
        qs = {q for q, x in zip(QUADS, m) if x}
        if sh.max() >= LINE_POINT:
            out.append("point")
        elif not qs:
            out.append("block")
        elif qs == {"--", "++"}:
            out.append("butterfly-inversion")
        elif qs == {"-+", "+-"}:
            out.append("butterfly-deldup")
        else:
            out.append("junction")
    return out


def detect_junctions(ui, uj, uc, cov_far5, n, bs, hot_params):
    """§9.3 steps 1-9 on the 5 kb upper triangle and the corrected covFar.

    Returns the model and one record per corner after NMS (unranked, untruncated; see junction_records).
    """
    tm, t = {}, time.time()

    def lap(name):
        nonlocal t
        tm[name] = round(time.time() - t, 2)
        t = time.time()

    jm = JunctionModel(ui, uj, uc, cov_far5, n, bs, hot_params)
    lap("model")
    I, J = jm.seed_boxes()
    lap("seeds")
    V = np.unique(jm.proxy_vertices(I, J), axis=0) if I.size else np.zeros((0, 2), np.int64)
    lap("proxy")
    sc = jm.score(V[:, 0], V[:, 1]) if len(V) else np.zeros(0)
    keep = sc >= TAU_SCREEN
    lap("screen")
    ca2, cb2, cl = jm.climb(V[keep, 0].copy(), V[keep, 1].copy(), sc[keep].copy())
    lap("climb")
    k = nms(ca2, cb2, cl)
    a2, b2, sc_p, unc = jm.polish(ca2[k], cb2[k])
    k = nms(a2, b2, sc_p)
    a2, b2, unc = a2[k], b2[k], unc[k]
    lap("polish")

    st = jm.quad_stats(a2, b2)
    score, mask = jm.choose(st)
    share = jm.line_share(a2, b2)
    ext = np.zeros((a2.size, 4), np.int64)
    n_merged = np.zeros(a2.size, np.int64)
    for k in range(a2.size):
        mem = np.maximum(np.abs(ca2 - a2[k]), np.abs(cb2 - b2[k])) <= NMS_HALF
        n_merged[k] = mem.sum()
        ext[k] = (ca2[mem].min(), ca2[mem].max(), cb2[mem].min(), cb2[mem].max()) if mem.any() else \
            (a2[k], a2[k], b2[k], b2[k])
    lap("final")
    jm.timing = tm
    jm.counts = dict(seeds=int(I.size), starts=len(V), screened=int(keep.sum()), corners=int(a2.size))
    return jm, dict(a2=a2, b2=b2, score=score, mask=mask, share=share, unc=unc, ext=ext, nMerged=n_merged,
                    type=classify(mask, share), **st)


def link_hotspots(A, B, hotspots, bs):
    """§9.3.11: ranks of hotspots whose bbox expanded by L bins contains (A, B), ascending, at most 20."""
    pad = JUN_L * bs
    return [h["rank"] for h in sorted(hotspots, key=lambda h: h["rank"])
            if h["y0"] - pad <= A < h["y1"] + pad and h["x0"] - pad <= B < h["x1"] + pad][:20]


def exp_fp(score, n_boxes):
    """§9.4: expected false corners per chromosome scoring >= score."""
    p = JUN_NULL["pU"] * np.exp(-(np.asarray(score, float) - JUN_NULL["u"]) / JUN_NULL["beta"])
    return JUN_NULL["basin"] * n_boxes * p


def junction_records(jm, cr, hotspots):
    """§9.3.10 + §9.6: corners with score >= TAU, ranked by score, as junctions document records."""
    bs = jm.bs
    keep = np.flatnonzero(cr["score"] >= TAU)
    keep = keep[np.argsort(-cr["score"][keep], kind="stable")]
    ef = exp_fp(cr["score"][keep], jm.eligible_boxes())
    q = np.minimum.accumulate((ef / np.arange(1, keep.size + 1))[::-1])[::-1] if keep.size else ef
    out = []
    for rank, (k, e, qq) in enumerate(zip(keep, ef, q), 1):
        A, B = int(cr["a2"][k] * bs // 2), int(cr["b2"][k] * bs // 2)
        quads = {}
        for jq, qd in enumerate(QUADS):
            ex = cr["exp"][k, jq]
            quads[qd] = dict(llr=sig4(cr["llr"][k, jq]), core=sig4(cr["core"][k, jq]), s=sig4(cr["s"][k, jq]),
                             obs=int(cr["obs"][k, jq]), exp=sig4(ex),
                             oe=sig4(cr["obs"][k, jq] / ex) if ex > 0 else None,
                             inner=sig4(cr["inner"][k, jq]), detected=bool(cr["mask"][k, jq]))
        det = [qd for jq, qd in enumerate(QUADS) if cr["mask"][k, jq]]
        typ = cr["type"][k]
        event = {"butterfly-inversion": "inversion", "butterfly-deldup": "deletion+duplication",
                 "junction": SINGLE_EVENT[det[0]] if det else None}.get(typ)
        y0, y1, x0, x1 = (int(v * bs // 2) for v in cr["ext"][k])
        sh = cr["share"][k]
        out.append(dict(rank=rank, A=A, B=B, type=typ, event=event, quads=det, score=sig4(cr["score"][k]),
                        expFP=sig4(max(e, FP_FLOOR)), q=sig4(max(qq, FP_FLOOR)),
                        vertexUncBp=int(cr["unc"][k] * bs // 2), extent=[y0, y1 + bs, x0, x1 + bs],
                        nMerged=int(cr["nMerged"][k]), lineShare=sig4(sh.max()),
                        lineAxis="row" if sh[0] >= sh[1] else "col", quadrants=quads,
                        hotspots=link_hotspots(A, B, hotspots, bs)))
    return out


def junctions_doc(jm, cr, hotspots, name, chrom_length):
    """§9.6 junctions document (embedded in the .bcm trailer)."""
    js = junction_records(jm, cr, hotspots)
    params = dict(resolution=jm.bs, L=JUN_L, dMinBp=jm.dmin * jm.bs, phiNear=sig4(jm.phi_near),
                  phiFar=sig4(jm.phi_far), template=[sig4(x) for x in jm.t], tau=TAU, tauScreen=TAU_SCREEN,
                  tau2=TAU2, coreMin=CORE_MIN, coreFrac=CORE_FRAC, innerD=INNER_D, innerSup=INNER_SUP,
                  lineSharePoint=LINE_POINT, nmsBp=NMS_HALF * jm.bs // 2, seedBox=JUN_G, seedObs=SEED_OBS,
                  seedOE=SEED_OE, nBoxes=jm.eligible_boxes(), null=JUN_NULL)
    counts = {t: sum(j["type"] == t for j in js) for t in JUN_TYPES}
    return dict(version=1, chrom=name, binSize=jm.bs, nBins=jm.n, chromLength=chrom_length, params=params,
                counts=counts, junctions=js)


def print_junctions(js, top=15):
    for j in js[:top]:
        print(f"  #{j['rank']:<3d} {j['A'] / 1e6:8.4f} {j['B'] / 1e6:8.4f} {j['type']:<19s} {'/'.join(j['quads']):<5s} "
              f"score {j['score']:>7} expFP {j['expFP']:.2g} q {j['q']:.2g} unc {j['vertexUncBp']} "
              f"line {j['lineShare']:.2f} hot {j['hotspots'][:4]}", flush=True)


def junction_null(jm, hotspots, n_sample, seed=0):
    """§9.4 optional re-calibration: score random eligible boxes like seeds; print the tail fit and expFP = 1."""
    t = time.time()
    n, G, bs = jm.n, JUN_G, jm.bs
    nb = n // G
    rng = np.random.default_rng(seed)
    I, J = rng.integers(0, nb, 4 * n_sample), rng.integers(0, nb, 4 * n_sample)
    ok = ((J - I) * G >= jm.dmin) & jm.valid[np.minimum(I * G, n - 1)] & jm.valid[np.minimum(J * G, n - 1)]
    I, J = I[ok][:n_sample], J[ok][:n_sample]
    # set aside boxes whose proxy window [x - R - L - 1, + win) overlaps a hotspot bbox padded by L + 8 bins
    win, back, pad = G + 2 * JUN_R + 2 * JUN_L + 2, JUN_R + JUN_L + 1, JUN_L + 8
    r0, c0 = I * G - back, J * G - back
    near = np.zeros(I.size, bool)
    for h in hotspots:
        near |= ((r0 < ceil_div(h["y1"], bs) + pad) & (r0 + win > h["y0"] // bs - pad)
                 & (c0 < ceil_div(h["x1"], bs) + pad) & (c0 + win > h["x0"] // bs - pad))
    V = jm.proxy_vertices(I, J)
    sc = jm.score(V[:, 0], V[:, 1])
    k = np.flatnonzero(sc >= TAU_SCREEN)
    if k.size:
        a2, b2, _ = jm.climb(V[k, 0].copy(), V[k, 1].copy(), sc[k].copy())
        sc[k] = jm.polish(a2, b2)[2]
    x, u = sc[~near], JUN_NULL["u"]
    n_boxes = jm.eligible_boxes()
    log(f"junction null: {I.size} boxes ({int(near.sum())} near hotspots set aside), {k.size} screened "
        f"({time.time() - t:.1f} s); quantiles 50/99/99.9/99.99% "
        + " ".join(f"{v:.1f}" for v in np.quantile(x, [.5, .99, .999, .9999])) + f" max {x.max():.1f}")
    e = x[x >= u] - u
    if e.size == 0:
        log(f"junction null: no exceedances over u = {u}")
        return
    beta, p_u = float(e.mean()), e.size / x.size
    tau1 = u + beta * math.log(JUN_NULL["basin"] * n_boxes * p_u)
    log(f"junction null: tail over u = {u}: n {e.size} beta {beta:.1f} pU {p_u:.3g}; nBoxes {n_boxes}, "
        f"basin {JUN_NULL['basin']} -> expFP = 1 at score {tau1:.0f} (stored constants: beta {JUN_NULL['beta']}, "
        f"pU {JUN_NULL['pU']})")

def run_pipeline(i, j, c, n, name, out, bin_size, far_bp, test_min_sep, chrom_length, junction_null_n=0):
    """§3.1 + §3.2 + §8.1 + §9 + §10.1: hotspots first (they correct covFar), then the contacts encoding, then the
    junction detector (run on the in-memory upper triangle and the rounded hotspot params, i.e. exactly the inputs
    --verify re-reads). Writes OUT/<name>.bcm and returns its path."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    t = time.time()
    check_symmetric(i, j, c, n, bin_size)
    params, hotspots, stripes, lowcov, cov_far = detect_hotspots(i, j, c, n, bin_size, far_bp, test_min_sep,
                                                                 chrom_length)
    log(f"hotspots done ({time.time() - t:.1f} s)")

    t = time.time()
    contacts, info = encode(i, j, c, n, bin_size, far_bp, cov_far)
    log(f"encode: nBins {n} nnz {info['nnz']} totalCount {info['totalCount']} diagCount {info['diagCount']} "
        f"nOverflow {info['nOverflow']} colBytes {info['colBytes']} contacts raw {len(contacts)} B "
        f"({time.time() - t:.1f} s)")

    doc = dict(version=1, chrom=name, binSize=bin_size, nBins=n, chromLength=chrom_length,
               nnz=info["nnz"], totalCount=info["totalCount"], params=params, hotspots=hotspots,
               stripes=stripes, lowcov=lowcov)

    t = time.time()
    up = i <= j
    jm, cr = detect_junctions(i[up], j[up], c[up], cov_far, n, bin_size, params)
    jdoc = junctions_doc(jm, cr, hotspots, name, chrom_length)
    log(f"junctions: phiNear {jm.phi_near:.4g} phiFar {jm.phi_far:.4g} {jm.counts} nBoxes "
        f"{jdoc['params']['nBoxes']} called {len(jdoc['junctions'])} {jdoc['counts']}")
    log(f"junctions timing (s): {jm.timing} ({time.time() - t:.1f} s)")
    print_junctions(jdoc["junctions"])
    if junction_null_n:
        junction_null(jm, hotspots, junction_null_n)

    meta = bundle_json(name, chrom_length, doc, jdoc)
    payload = bundle_payload(contacts, info["colBytes"], meta)
    gz = gzip.compress(payload, compresslevel=9, mtime=0)
    path = out / f"{name}.bcm"
    path.write_bytes(gz)
    log(f"wrote {path}: payload {len(payload)} B (contacts {len(contacts)} B, json {len(meta)} B), gz {len(gz)} B")
    return path


# ------------------------------------------------------------- §7.1 verify

class Checks:
    """Collect pass/fail results for --verify."""

    def __init__(self):
        self.failed = 0

    def __call__(self, ok, msg):
        print(f"  {'ok  ' if ok else 'FAIL'} {msg}", flush=True)
        self.failed += not ok
        return bool(ok)


def varint_decode(buf, count):
    """Decode `count` unsigned LEB128 varints from a uint8 array (vectorised)."""
    end = (buf & 0x80) == 0
    if int(end.sum()) != count or (buf.size and not end[-1]):
        raise ValueError(f"cols stream holds {int(end.sum())} varints, expected {count}")
    if count == 0:
        return np.zeros(0, dtype=np.int64)
    vid = np.concatenate(([0], np.cumsum(end)[:-1]))
    starts = np.concatenate(([0], np.flatnonzero(end)[:-1] + 1))
    pos = np.arange(buf.size) - starts[vid]
    vals = (buf & 0x7F).astype(np.uint64) << (7 * pos).astype(np.uint64)
    return np.bincount(vid, weights=vals.astype(np.float64), minlength=count).astype(np.int64)


def compact_size(doc):
    """Byte length of a document serialised as in the bundle trailer."""
    return len(json.dumps(doc, separators=(",", ":")).encode("utf-8"))


class BundleError(Exception):
    """A .bcm file that cannot be decoded; the message is shown as a FAIL line."""


class PayloadReader:
    """Bounded reads from a gzip stream (or a raw payload), so a forged header cannot force a huge allocation."""

    def __init__(self, raw):
        self.gz = raw[:2] == b"\x1f\x8b"
        self.raw, self.pos = raw, 0
        self.d = zlib.decompressobj(wbits=31) if self.gz else None
        self.tail = raw

    def read(self, k):
        """Exactly k payload bytes, else BundleError (truncated or corrupt)."""
        if not self.gz:
            out = self.raw[self.pos:self.pos + k]
        else:
            parts, need = [], k
            try:
                while need > 0 and not self.d.eof:
                    chunk = self.d.decompress(self.tail, need)
                    self.tail = self.d.unconsumed_tail
                    parts.append(chunk)
                    need -= len(chunk)
                    if not chunk and not self.tail:
                        break
            except zlib.error as e:
                raise BundleError(f"gzip stream is corrupt: {e}") from None
            out = b"".join(parts)
        self.pos += len(out)
        if len(out) != k:
            raise BundleError(f"truncated or corrupt: payload ends at {self.pos} B, expected at least "
                              f"{self.pos - len(out) + k} B")
        return out

    def end_error(self):
        """None if the payload ends exactly here (gzip: one complete member, nothing after it), else why not."""
        if not self.gz:
            return None if self.pos == len(self.raw) else f"{len(self.raw) - self.pos} B after the expected end"
        try:
            extra = self.d.decompress(self.tail, 1)
        except zlib.error as e:
            return f"gzip stream is corrupt: {e}"
        if extra:
            return "payload continues after the expected end"
        if not self.d.eof:
            return "gzip stream is incomplete (truncated)"
        if self.d.unused_data:
            return f"{len(self.d.unused_data)} B of trailing data after the gzip stream"
        return None


def decode_bundle(raw):
    """Reference decoder for §10.1 (sections as §2.2). Returns header dict, level-0 arrays and the trailer JSON.

    Raises BundleError with a readable message for anything that is not a well-formed bundle. The header is
    checked against MAX_NBINS / MAX_NNZ and the section sizes are known before any section is read.
    """
    rd = PayloadReader(raw)
    hdr = dict(gzip=rd.gz, magic=raw[:4] if not rd.gz else None, version=None)
    head = rd.read(HEADER.size)
    magic, ver, bs, n, nnz, nov, colb, minsep, total, diag = HEADER.unpack(head)
    hdr.update(magic=magic, version=ver, binSize=bs, nBins=n, nnz=nnz, nOverflow=nov, colBytes=colb,
               minSepBp=minsep, totalCount=total, diagCount=diag)
    if magic == b"BCM1":
        raise BundleError("BCM1 file (old contacts.bin format); re-run preprocess.py")
    if magic != MAGIC or ver != VERSION:
        raise BundleError(f"bad magic {magic!r} / version {ver}")
    limits = [(bs > 0, f"binSize {bs} > 0"), (n <= MAX_NBINS, f"nBins {n} <= {MAX_NBINS}"),
              (nnz <= MAX_NNZ, f"nnz {nnz} <= {MAX_NNZ}"),
              (nnz <= n * (n + 1) // 2, f"nnz {nnz} <= nBins(nBins+1)/2 = {n * (n + 1) // 2}"),
              (nov <= nnz, f"nOverflow {nov} <= nnz"),
              (nnz <= colb <= 5 * nnz, f"colBytes {colb} within [nnz, 5 nnz] (1-5 byte varints)")]
    bad = [msg for ok, msg in limits if not ok]
    if bad:
        raise BundleError("header out of range: " + "; ".join(bad))
    pad, cpad = (-nnz) % 4, (-colb) % 4
    sections = 4 * (n + 1) + 8 * n + nnz + pad + 4 * nov + colb + cpad
    body = rd.read(sections + 4)
    json_len = U32.unpack_from(body, sections)[0]
    hdr["jsonLen"], hdr["expectLen"] = json_len, 64 + sections + 4 + json_len
    if json_len > MAX_TRAILER:
        raise BundleError(f"jsonLen {json_len} > {MAX_TRAILER}")
    text = rd.read(json_len)
    err = rd.end_error()
    if err:
        raise BundleError(f"{err} (expected payload {hdr['expectLen']} B)")
    hdr["payloadLen"] = hdr["expectLen"]
    try:
        meta = json.loads(text.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as e:
        raise BundleError(f"trailer is not UTF-8 JSON: {e}") from None
    if not isinstance(meta, dict):
        raise BundleError("trailer JSON is not an object")

    o = 0
    row_ptr = np.frombuffer(body, "<u4", n + 1, o).astype(np.int64); o += 4 * (n + 1)
    cov_all = np.frombuffer(body, "<u4", n, o).astype(np.int64); o += 4 * n
    cov_far = np.frombuffer(body, "<u4", n, o).astype(np.int64); o += 4 * n
    counts = np.frombuffer(body, np.uint8, nnz, o).astype(np.int64); o += nnz
    hdr["padZero"] = body[o:o + pad] == bytes(pad); o += pad
    overflow = np.frombuffer(body, "<u4", nov, o).astype(np.int64); o += 4 * nov
    cols = np.frombuffer(body, np.uint8, colb, o); o += colb
    hdr["padZero"] &= body[o:o + cpad] == bytes(cpad)
    if row_ptr[0] != 0 or row_ptr[-1] != nnz or np.any(np.diff(row_ptr) < 0):
        raise BundleError(f"rowPtr is not a CSR pointer array (first {row_ptr[0]}, last {row_ptr[-1]}, nnz {nnz})")
    esc = counts == 255
    hdr["escapes"] = int(esc.sum())
    if hdr["escapes"] != nov:
        raise BundleError(f"{hdr['escapes']} escaped counts but nOverflow {nov}")
    counts[esc] = overflow
    end = (cols & 0x80) == 0
    if nnz and int(np.diff(np.concatenate(([-1], np.flatnonzero(end)))).max(initial=1)) > 5:
        raise BundleError("cols holds a varint longer than 5 bytes")
    try:
        delta = varint_decode(cols, nnz)
    except ValueError as e:
        raise BundleError(str(e)) from None
    row_len = np.diff(row_ptr)
    rows = np.repeat(np.arange(n), row_len)
    cs = np.cumsum(delta)
    starts = row_ptr[:-1][row_len > 0]
    base = np.repeat(cs[starts] - delta[starts], row_len[row_len > 0])
    cols_j = rows + cs - base
    if nnz and (cols_j.min() < 0 or cols_j.max() >= n):
        raise BundleError(f"decoded columns out of range [0, {n})")
    return hdr, dict(row_ptr=row_ptr, i=rows, j=cols_j, c=counts, delta=delta, first=starts,
                     cov_all=cov_all, cov_far=cov_far), meta


def verify(path, input_path):
    """--verify entry point: any unexpected error in a malformed bundle becomes a FAIL line and exit status 1."""
    try:
        return verify_bundle(path, input_path)
    except Exception as e:  # noqa: BLE001 - a verifier must not crash on hostile input
        print(f"  FAIL verify aborted: {type(e).__name__}: {e}", flush=True)
        print("verify: FAIL")
        return 1


def verify_bundle(path, input_path):
    """§7.1 + §8.1 + §9.7 + §10.1: re-decode a .bcm bundle and check it against the TSV and its own JSON."""
    path = Path(path)
    check = Checks()
    try:
        raw = path.read_bytes()
    except OSError as e:
        check(False, f"cannot read {path}: {e.strerror}")
        print("verify: FAIL")
        return 1
    print(f"verify {path}: {len(raw)} B")
    try:
        hdr, m, meta = decode_bundle(raw)
    except BundleError as e:
        check(False, f"decode: {e}")
        print("verify: FAIL")
        return 1
    check(hdr["gzip"], "file is gzip-compressed")
    check(True, f"magic {hdr['magic']!r} version {hdr['version']}, header within limits, payload length "
          f"{hdr['payloadLen']} == expected")
    check(hdr["padZero"], "padding after counts and after cols is zero")
    hot, jun = meta.get("hotspots"), meta.get("junctions")
    n, bs = hdr["nBins"], hdr["binSize"]
    print(f"  json trailer {hdr['jsonLen']} B: hotspots {compact_size(hot)} B, junctions "
          f"{compact_size(jun) if jun is not None else 'null'}")
    ok = (meta.get("format") == "bcm" and meta.get("version") == VERSION and isinstance(meta.get("chrom"), str)
          and bool(meta.get("chrom")) and has_type(meta.get("chromLength"), int) and isinstance(hot, dict)
          and (jun is None or isinstance(jun, dict)) and "createdBy" in meta)
    check(ok, f"trailer keys: format {meta.get('format')!r} version {meta.get('version')} "
          f"chrom {meta.get('chrom')!r} chromLength {meta.get('chromLength')}")
    cl = meta.get("chromLength")
    ok &= check(isinstance(cl, int) and cl > (n - 1) * bs and n == ceil_div(cl, bs),
                f"nBins {n} == ceil(chromLength {cl} / binSize {bs})")
    docs = [d for d in (hot, jun) if isinstance(d, dict)]
    ok &= check(all(d.get("chrom") == meta.get("chrom") and d.get("chromLength") == cl and d.get("binSize") == bs
                    and d.get("nBins") == n for d in docs),
                "chrom, chromLength, binSize, nBins agree between trailer, header and the embedded documents")
    if not ok:
        print(f"verify: FAIL ({check.failed} checks; trailer inconsistent, later checks skipped)")
        return 1
    print(f"  header: binSize {bs} nBins {n} nnz {hdr['nnz']} nOverflow {hdr['nOverflow']} colBytes "
          f"{hdr['colBytes']} minSepBp {hdr['minSepBp']} totalCount {hdr['totalCount']:.0f} "
          f"diagCount {hdr['diagCount']:.0f}")
    i, j, c = m["i"], m["j"], m["c"]
    check(hdr["minSepBp"] == FAR_BP, f"header minSepBp {hdr['minSepBp']} == {FAR_BP} (covFar definition)")
    check(hdr["escapes"] == hdr["nOverflow"], f"escape count {hdr['escapes']} == nOverflow")
    check(int((c >= 255).sum()) == hdr["nOverflow"], "nOverflow == #cells with count >= 255")
    check(bool(np.all(m["delta"][m["first"]] >= 0)) and
          bool(np.all(np.delete(m["delta"], m["first"]) >= 1)), "col deltas: first >= 0, later >= 1")
    check(bool(np.all(j >= i)) and bool(np.all(j < n)) and bool(np.all(c > 0)), "entries upper, in range, nonzero")
    off = i != j
    total, diag = int(2 * c[off].sum() + c[~off].sum()), int(c[~off].sum())
    check(total == hdr["totalCount"] and diag == hdr["diagCount"], f"header total/diag match decoded ({total}, {diag})")

    mi = np.concatenate((i, j[off]))
    mj = np.concatenate((j, i[off]))
    mc = np.concatenate((c, c[off]))
    far_bins = hdr["minSepBp"] // bs
    cov_all, cov_far_raw = coverages(mi, mj, mc, n, far_bins)
    check(np.array_equal(cov_all, m["cov_all"]), "covAll == decoded row sums")
    removed = cov_far_raw - m["cov_far"]
    check(bool(np.all(removed >= 0)), f"covFar <= raw far coverage from decoded (self-excluded: "
          f"{int((removed > 0).sum())} bins, {int(removed.sum())} counts)")

    tsv = None
    if input_path:
        tsv = load_tsv(input_path, bs)
    elif meta.get("chrom") == SYN_NAME:
        si, sj, sc, _ = synthetic_tsv()
        tsv = (si, sj, sc)
        print("  (comparing against regenerated synthetic TSV)")
    if tsv is not None:
        ti, tj, tc = tsv
        t_n = int(max(ti.max(), tj.max())) + 1
        tup = ti <= tj
        t_total, t_diag = int(tc.sum()), int(tc[ti == tj].sum())
        # --fai / --chrom-length may add empty bins after the last TSV bin
        check(t_n == n or (t_n < n and cl > (n - 1) * bs), f"nBins {n} == TSV {t_n} (or chromLength extends it)")
        check(int(tup.sum()) == hdr["nnz"], f"nnz {hdr['nnz']} == TSV {int(tup.sum())}")
        check(t_total == hdr["totalCount"], f"totalCount == TSV {t_total}")
        check(t_diag == hdr["diagCount"], f"diagCount == TSV {t_diag}")
        check(int((tc[tup] >= 255).sum()) == hdr["nOverflow"], f"nOverflow == TSV {int((tc[tup] >= 255).sum())}")
        ka, kb = mi * n + mj, ti.astype(np.int64) * n + tj
        oa, ob = np.argsort(ka), np.argsort(kb)
        same = ka.size == kb.size and np.array_equal(ka[oa], kb[ob]) and np.array_equal(mc[oa], tc[ob])
        check(same, f"decoded mirrored matrix == all {kb.size} TSV lines")
        t_all, t_far = coverages(ti, tj, tc, n, far_bins)
        check(np.array_equal(t_all, m["cov_all"]), "covAll == TSV bincount")
        check(bool(np.all(m["cov_far"] <= t_far)), "covFar <= TSV raw far coverage")
        del tsv, ti, tj, tc, ka, kb, oa, ob

    for L in range(N_LEVELS):
        I, J = i >> L, j >> L
        nL = ceil_div(n, 1 << L)
        w = np.where((I == J) & off, 2 * c, c)
        if nL <= 4096:
            C = np.bincount(I * nL + J, weights=w, minlength=nL * nL).reshape(nL, nL)
            dsum, osum = float(np.trace(C)), float(np.triu(C, 1).sum())
        else:
            dsum, osum = float(w[I == J].sum()), float(w[I < J].sum())
        cov_l = np.bincount(np.arange(n) >> L, weights=m["cov_all"], minlength=nL)
        ok = 2 * osum + dsum == hdr["totalCount"] and cov_l.sum() == hdr["totalCount"]
        check(ok, f"level {L} ({bs << L} bp, n {nL}): 2*upper + diag == totalCount == sum covAll")
    check(len(raw) <= MAX_BIN_SIZE, f"{path.name} {len(raw)} B <= {MAX_BIN_SIZE}")

    verify_hotspots(hot, hdr, check, (mi, mj, mc), m["cov_far"], (i, j, c))
    verify_junctions(jun, hot, hdr, check, (i, j, c), m["cov_far"])
    print(f"verify: {'PASS' if check.failed == 0 else f'FAIL ({check.failed} checks)'}")
    return 1 if check.failed else 0


def verify_hotspots(hot, hdr, check, full, cov_far, upper):
    """§7.1 hotspot checks, §8.1 fields and params, plus the fixture expectations of §3.3."""
    p = hot["params"]
    res, minsep = p["resolution"], p["minSepBp"]
    for k in ("binSize", "nBins", "nnz", "totalCount"):
        check(hot[k] == hdr[k], f"hotspots {k} {hot[k]} == header")
    check(minsep == TEST_MIN_SEP and p["lowcovFactor"] == LOWCOV_FACTOR,
          f"params.minSepBp {minsep}, lowcovFactor {p['lowcovFactor']}")
    pbl = p.get("phiByLevel", [])
    # ebalFar50k is 0 only when nothing is testable (no valid far cells), e.g. an empty or tiny chromosome
    check(len(pbl) == N_LEVELS and all(x >= 0 for x in pbl) and (p.get("ebalFar50k", 0) > 0 or p["nTests"] == 0)
          and p.get("phiNear", -1) >= 0 and p.get("phiFar", -1) >= 0,
          f"params phiNear {p.get('phiNear')} phiFar {p.get('phiFar')} ebalFar50k {p.get('ebalFar50k')} "
          f"phiByLevel {pbl}")
    re_pbl = phi_by_level(*upper, cov_far, hdr["nBins"], hdr["binSize"], p["ebalFar50k"])
    check(len(pbl) == N_LEVELS and np.allclose(re_pbl, pbl, rtol=1e-3, atol=1e-6),
          f"phiByLevel recomputed from decoded data == params ({[sig4(x) for x in re_pbl]})")

    hs = hot["hotspots"]
    size = compact_size(hot)
    check(size <= MAX_JSON_SIZE, f"hotspots document {size} B <= {MAX_JSON_SIZE}")
    check([h["rank"] for h in hs] == list(range(1, len(hs) + 1)), f"ranks contiguous 1..{len(hs)}")
    check(all(a["excess"] >= b["excess"] for a, b in zip(hs, hs[1:])), "sorted by excess descending")
    # §7.1's bbox rule (x0 - y1 >= minSepBp - resolution) cannot hold once near-band cells (300 kb) merge
    # into components; every tested cell, and so every peak, is still >= minSepBp off the diagonal.
    bad_peak = [h["rank"] for h in hs if h["peakX"] - h["peakY"] < minsep]
    check(not bad_peak, f"peakX - peakY >= minSepBp = {minsep} for all (violations: {bad_peak[:10]})")
    bad_sep = [(h["rank"], h["x0"] - h["y1"]) for h in hs if h["x0"] - h["y1"] < minsep - res]
    print(f"  note {len(bad_sep)} merged boxes have x0 - y1 < {minsep - res} (rank, gap bp): {bad_sep[:10]}")
    bad_y = [h["rank"] for h in hs if h["y1"] > h["x0"]]
    check(not bad_y, f"y1 <= x0 for all (violations: {bad_y[:10]})")

    def hits(h, s):
        return (h["y0"] < s["end"] and s["start"] < h["y1"]) or (h["x0"] < s["end"] and s["start"] < h["x1"])
    bad_stripe = [h["rank"] for h in hs if any(hits(h, s) for s in hot["stripes"])]
    check(not bad_stripe, f"no hotspot intersects a stripe row/col (violations: {bad_stripe[:10]})")

    fields = ("nz5k", "maxRowShare5k", "maxColShare5k", "hubDegree", "shape")
    check(all(all(k in h for k in fields) for h in hs), "every hotspot has " + ", ".join(fields))
    re = [{k: h[k] for k in ("y0", "y1", "x0", "x1")} for h in hs]
    shape_fields(re, *full, hot["binSize"])
    hub = hub_degrees(hs)
    bad = [h["rank"] for h, r, k in zip(hs, re, hub)
           if any(h[f] != r[f] for f in fields[:3] + ("shape",)) or h["hubDegree"] != k]
    check(not bad, f"5 kb shape fields and hubDegree recomputed from decoded data match (mismatch: {bad[:10]})")

    if hot["chrom"] == SYN_NAME:
        bs = hot["binSize"]

        def near(h, ev):
            y, x = ev[0] * bs, ev[1] * bs
            return (h["y0"] - 100_000 <= y < h["y1"] + 100_000) and (h["x0"] - 100_000 <= x < h["x1"] + 100_000)
        top = hs[:2]
        ok = len(top) == 2 and any(near(top[0], a) and near(top[1], b)
                                   for a, b in (SYN_EVENTS, SYN_EVENTS[::-1]))
        check(ok, "fixture: hotspots 1-2 match injected events (3.0 x 4.5 Mb, 6.0 x 8.5 Mb) within 100 kb")
        sy = SYN_STRIPE[0] * bs
        check(any(s["start"] <= sy < s["end"] for s in hot["stripes"]), "fixture: stripe at 7.5 Mb")
        check(any(a <= 500_000 and b >= 700_000 for a, b in hot["lowcov"]), "fixture: lowcov covers 0.5-0.7 Mb")

    print(f"  {len(hs)} hotspots ({sum(h['shape'] == 'point' for h in hs)} point), {len(hot['stripes'])} stripes, "
          f"{len(hot['lowcov'])} lowcov runs; nTests {p['nTests']} phiNear {p['phiNear']} phiFar {p['phiFar']} "
          f"ebalFar50k {p['ebalFar50k']}")
    print("  rank  y (Mb)          x (Mb)            obs      exp    excess  peakOE  -log10p cells  nz5k "
          "rowSh  colSh  shape     hub")
    for h in hs[:15]:
        print(f"  {h['rank']:>4}  {h['y0'] / 1e6:6.2f}-{h['y1'] / 1e6:6.2f}  {h['x0'] / 1e6:6.2f}-{h['x1'] / 1e6:6.2f}"
              f"  {h['obs']:>6} {h['exp']:>8} {h['excess']:>9} {str(h['peakOE']):>7} {h['peakNegLog10P']:>8}"
              f" {h['nCells']:>5} {h['nz5k']:>5} {h['maxRowShare5k']:>6} {h['maxColShare5k']:>6}"
              f"  {h['shape']:<9} {h['hubDegree']:>3}")
    print("  stripes: " + ", ".join(f"{s['start'] / 1e6:.2f}-{s['end'] / 1e6:.2f} (x{s['maxFactor']})"
                                    for s in hot["stripes"]))


JUN_KEYS = dict(rank=int, A=int, B=int, type=str, event=(str, type(None)), quads=list, score=float, expFP=float,
                q=float, vertexUncBp=int, extent=list, nMerged=int, lineShare=float, lineAxis=str, quadrants=dict,
                hotspots=list)
QUAD_KEYS = dict(llr=float, core=float, s=float, obs=int, exp=float, oe=(float, type(None)), inner=float,
                 detected=bool)
JUN_PARAM_KEYS = ("resolution", "L", "dMinBp", "phiNear", "phiFar", "template", "tau", "tauScreen", "tau2",
                  "coreMin", "coreFrac", "innerD", "innerSup", "lineSharePoint", "nmsBp", "seedBox", "seedObs",
                  "seedOE", "nBoxes", "null")
TYPE_QUADS = {"butterfly-inversion": ({"--", "++"}, "inversion"), "butterfly-deldup": ({"-+", "+-"},
                                                                                      "deletion+duplication")}


def has_type(x, t):
    """JSON type check; float accepts ints, int rejects bools."""
    ts = t if isinstance(t, tuple) else (t,)
    if isinstance(x, bool):
        return bool in ts
    return any(isinstance(x, (int, float) if tt is float else tt) for tt in ts)


def verify_junctions(doc, hot, hdr, check, upper, cov_far):
    """§9.7: junctions document schema, consistency, hotspot links, determinism and the fixture expectations."""
    if doc is None:
        print("  note: bundle has junctions = null, junction checks skipped")
        return
    js, p = doc.get("junctions", []), doc.get("params", {})
    size = compact_size(doc)
    print(f"  junctions document {size} B, {len(js)} junctions")
    check(size <= MAX_JSON_SIZE, f"junctions document {size} B <= {MAX_JSON_SIZE}")
    top = all(k in doc for k in ("version", "chrom", "binSize", "nBins", "chromLength", "params", "counts",
                                 "junctions"))
    check(top and all(k in p for k in JUN_PARAM_KEYS) and len(p.get("template", [])) == p.get("L"),
          "junctions top-level keys and params present")
    check(all(doc.get(k) == hot.get(k) for k in ("chrom", "binSize", "nBins", "chromLength")),
          "chrom, binSize, nBins, chromLength match the hotspots document")

    bad = []
    for j in js:
        ok = all(k in j and has_type(j[k], t) for k, t in JUN_KEYS.items())
        ok = ok and set(j["quadrants"]) == set(QUADS) and all(
            all(k in qd and has_type(qd[k], t) for k, t in QUAD_KEYS.items()) for qd in j["quadrants"].values())
        ok = ok and len(j["extent"]) == 4 and all(has_type(x, int) for x in j["extent"])
        ok = ok and j["type"] in JUN_TYPES and j["lineAxis"] in ("row", "col") and len(j["hotspots"]) <= 20
        if not ok:
            bad.append(j.get("rank"))
    check(not bad, f"every junction has the §9.6 keys and types (bad ranks: {bad[:10]})")
    if bad:
        return
    check([j["rank"] for j in js] == list(range(1, len(js) + 1)), f"ranks contiguous 1..{len(js)}")
    check(all(a["score"] >= b["score"] for a, b in zip(js, js[1:])) and all(j["score"] >= p["tau"] for j in js),
          f"scores sorted descending and >= tau {p['tau']}")
    check(all(j["B"] > j["A"] and j["B"] - j["A"] >= p["dMinBp"] for j in js), f"A < B and B - A >= {p['dMinBp']}")
    check(doc["counts"] == {t: sum(j["type"] == t for j in js) for t in JUN_TYPES}, f"counts {doc['counts']}")

    slack = 1 - 1e-3  # values in the file are rounded to 4 significant digits

    def supported(qd):
        return (qd["llr"] >= p["tau2"] * slack and qd["core"] >= max(p["coreMin"], p["coreFrac"] * qd["llr"]) * slack
                and qd["inner"] >= p["innerSup"] * slack)

    bad = []
    for j in js:
        qs = j["quads"]
        det = [q for q in QUADS if j["quadrants"][q]["detected"]]
        ok = sorted(qs) == sorted(det) and len(set(qs)) == len(qs) and all(supported(j["quadrants"][q]) for q in qs)
        if j["type"] in TYPE_QUADS:
            ok &= set(qs) == TYPE_QUADS[j["type"]][0] and len(qs) == 2 and j["event"] == TYPE_QUADS[j["type"]][1]
        elif j["type"] == "junction":
            ok &= len(qs) == 1 and j["event"] == SINGLE_EVENT[qs[0]]
        else:
            ok &= j["event"] is None and (j["type"] == "point" or not qs)
        if not ok:
            bad.append(j["rank"])
    check(not bad, f"type, event, quads, detected flags and support rule consistent (bad ranks: {bad[:10]})")

    ranks = {h["rank"] for h in hot["hotspots"]}
    bad = [j["rank"] for j in js if not set(j["hotspots"]) <= ranks
           or j["hotspots"] != link_hotspots(j["A"], j["B"], hot["hotspots"], doc["binSize"])]
    check(not bad, f"hotspot links are valid ranks and recompute identically (bad: {bad[:10]})")

    t = time.time()
    jm, cr = detect_junctions(*upper, cov_far, hdr["nBins"], hdr["binSize"], hot["params"])
    re = junction_records(jm, cr, hot["hotspots"])
    same = len(re) == len(js) and all((a["A"], a["B"], a["type"], a["quads"]) == (b["A"], b["B"], b["type"], b["quads"])
                                      for a, b in zip(js, re))
    rel = max((abs(a["score"] - b["score"]) / max(abs(a["score"]), 1e-12) for a, b in zip(js, re)), default=0.0)
    check(same and rel <= 1e-3, f"determinism: re-run from decoded data gives {len(re)} junctions, identical "
          f"(A, B, type, quads), max relative score difference {rel:.2g} <= 1e-3 ({time.time() - t:.1f} s)")

    if doc["chrom"] == SYN_NAME:
        bfly = [j for j in js if j["type"] == "butterfly-inversion" and sorted(j["quads"]) == ["++", "--"]
                and abs(j["A"] - 3_000_000) <= 10_000 and abs(j["B"] - 4_500_000) <= 10_000]
        # injected +- cells (1200+u, 1700-v): boundaries are exactly row 1200, column 1701
        jun = [j for j in js if j["type"] == "junction" and j["quads"] == ["+-"] and j["event"] == "tandem-dup"
               and j["A"] == 6_000_000 and j["B"] == 8_505_000]
        check(len(js) == 2 and len(bfly) == 1 and len(jun) == 1,
              "fixture: exactly 2 junctions, butterfly-inversion --/++ within 10 kb of 3.0/4.5 Mb and +- tandem-dup "
              "exactly at 6.0/8.505 Mb")
        check(bool(bfly and jun) and bfly[0]["hotspots"] == [1] and jun[0]["hotspots"] == [2],
              "fixture: butterfly links hotspot rank 1, junction links rank 2")
    else:
        print(f"  junction counts: {doc['counts']}")
        print_junctions(js)


# ---------------------------------------------------------------------- CLI

def warn(msg):
    print(f"warning: {msg}", file=sys.stderr, flush=True)


def chrom_name(path, fai=None):
    """§10.2 chromosome name from an input path.

    With a .fai index: the longest sequence name that is a prefix of the basename and is followed by '_', '.' or
    the end (NC_000001.11_whole.stat -> NC_000001.11). Otherwise, or if no fai name matches: the basename up to
    the first '_' or '.' (chr1_whole_5kb_bin.stat -> chr1).
    """
    base = Path(path).name
    if fai:
        hits = [k for k in fai if k and base.startswith(k) and (len(base) == len(k) or base[len(k)] in "_.")]
        if hits:
            return max(hits, key=len)
    return re.split(r"[_.]", base, maxsplit=1)[0]


def check_name(name, ap):
    if not name or name.startswith(".") or "/" in name or "\\" in name or "\0" in name:
        ap.error(f"chromosome name {name!r} cannot be used as a file name (give --name once per --input)")


def read_fai(path):
    """samtools faidx index: name<TAB>length<TAB>... per line. Returns {name: length}."""
    out = {}
    for k, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 2 or not f[1].strip().isdigit():
            sys.exit(f"{path}:{k}: not a .fai line (expected name<TAB>length): {line[:80]!r}")
        out[f[0]] = int(f[1])
    return out


def chrom_extent(name, n_data, bs, fai, chrom_length):
    """§10.2: (nBins, chromLength) from --chrom-length, else the .fai entry, else the data (nBins * binSize).

    A given length must cover every data bin: ceil(length / bs) >= n_data, i.e. length > (n_data - 1) * bs.
    """
    if chrom_length is not None:
        length, src = chrom_length, "--chrom-length"
    elif fai is not None and name in fai:
        length, src = fai[name], "--fai"
    else:
        if fai is not None:
            warn(f"{name!r} not found in the --fai index; using chromLength = nBins * binSize = {n_data * bs}")
        return n_data, n_data * bs
    if length <= (n_data - 1) * bs:
        sys.exit(f"{name}: chromLength {length} from {src} is shorter than the data, whose last bin starts at "
                 f"{(n_data - 1) * bs} (nBins {n_data} x binSize {bs})")
    n = ceil_div(length, bs)
    log(f"{name}: chromLength {length} from {src}; nBins {n} (data {n_data})")
    return n, length


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="Note: " + CALIBRATION_NOTE + ".")
    ap.add_argument("--input", nargs="+", metavar="STAT", help="one or more `a b count` files, one chromosome each")
    ap.add_argument("--name", action="append", metavar="NAME",
                    help="chromosome name; repeat once per --input, in order (default: from --fai or the file name)")
    ap.add_argument("--out", metavar="DIR", help="output directory; writes DIR/<name>.bcm per input")
    ap.add_argument("--fai", metavar="REF.fai", help="samtools faidx index giving the chromosome lengths")
    ap.add_argument("--bin", type=int, default=5000, help="bin size in bp of the input positions (default 5000)")
    ap.add_argument("--min-sep", type=int, default=FAR_BP, help="far-coverage / header minSepBp (bp)")
    ap.add_argument("--test-min-sep", type=int, default=TEST_MIN_SEP, help="hotspot testing separation (bp)")
    ap.add_argument("--chrom-length", type=int, help="chromosome length in bp (single --input; overrides --fai)")
    ap.add_argument("--synthetic", action="store_true", help="write the synthetic fixture DIR/synthetic.bcm")
    ap.add_argument("--verify", metavar="FILE.bcm", help="check a bundle (against --input STAT if given)")
    ap.add_argument("--junction-null", type=int, default=0, metavar="N",
                    help="re-run the random-box null that calibrates junction scores, with N boxes, and print the "
                         "score at which one false junction is expected")
    a = ap.parse_args()

    if a.verify:
        if a.input and len(a.input) > 1:
            ap.error("--verify takes at most one --input")
        if a.name:
            ap.error("--name cannot be used with --verify (the bundle carries its chromosome name)")
        sys.exit(verify(a.verify, a.input[0] if a.input else None))
    if not a.out or not (a.synthetic or a.input):
        ap.error("need --out with --input or --synthetic (or --verify FILE.bcm)")
    if a.synthetic and (a.input or a.name or a.fai or a.chrom_length is not None):
        ap.error("--synthetic cannot be combined with --input, --name, --fai or --chrom-length")
    if a.min_sep % a.bin or RESOLUTION % a.bin:
        ap.error("--min-sep and the 50 kb resolution must be multiples of --bin")
    t0 = time.time()
    if a.synthetic:
        i, j, c, n = synthetic_tsv()
        run_pipeline(i, j, c, n, SYN_NAME, a.out, SYN_BIN, a.min_sep, a.test_min_sep, n * SYN_BIN,
                     a.junction_null)
        log(f"done in {time.time() - t0:.1f} s")
        return

    if len(a.input) > 1 and a.chrom_length is not None:
        ap.error("--chrom-length is only allowed with a single --input (use --fai for several)")
    if a.name and len(a.name) != len(a.input):
        ap.error(f"--name given {len(a.name)} time(s) for {len(a.input)} input(s); give it once per --input, "
                 f"in order, or not at all")
    missing = [x for x in a.input if not Path(x).is_file()]
    if missing:
        ap.error(f"input file(s) not found: {missing}")
    fai = read_fai(a.fai) if a.fai else None
    names = a.name or [chrom_name(x, fai) for x in a.input]
    for name in names:
        check_name(name, ap)
    dup = sorted({x for x in names if names.count(x) > 1})
    if dup:
        ap.error(f"inputs map to the same chromosome name(s) {dup}, so each would overwrite DIR/<name>.bcm; "
                 f"pass --fai so names are matched against the index, or give --name once per --input")
    if not a.junction_null:
        warn(CALIBRATION_NOTE)

    written = []
    for path, name in zip(a.input, names):
        t = time.time()
        log(f"== {name} ({path})")
        i, j, c = load_tsv(path, a.bin)
        n_data = int(max(i.max(), j.max())) + 1
        log(f"load: {i.size} lines, nBins {n_data} ({time.time() - t:.1f} s)")
        n, length = chrom_extent(name, n_data, a.bin, fai, a.chrom_length)
        written.append(run_pipeline(i, j, c, n, name, a.out, a.bin, a.min_sep, a.test_min_sep, length,
                                    a.junction_null))
        del i, j, c
        log(f"{name} done in {time.time() - t:.1f} s")
    log(f"done in {time.time() - t0:.1f} s: " + ", ".join(str(x) for x in written))


if __name__ == "__main__":
    main()
