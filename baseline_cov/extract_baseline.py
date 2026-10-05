#!/usr/bin/env python3
"""Uncorrected window means and baseline means for the kept N400 trials.

Baseline as a regression covariate (Alday 2019, Psychophysiology 56:e13451)
instead of subtraction. The kept-trial set is the existing trials.csv, so the
artifact rule (on the baseline-corrected CP ptp) is unchanged. Each row is
checked against the stored dv: raw window - baseline must equal dv_cp / dv_fz.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import n400_common as C

OUT = C.OUT_ROOT / "baseline_cov"


def subject_rows(subject: str, onsets: dict, cp_ix: list[int], fz_ix: int) -> pd.DataFrame:
    tr = pd.read_csv(C.OUT_ROOT / subject / "trials.csv")
    tr = tr.loc[tr["keep"] == 1].copy()
    times = C.epoch_times()
    n_t = len(times)
    pre = int(round(-C.TMIN * C.SFREQ))
    base = times < -1e-8
    win = (times >= C.N400_WIN[0] - 1e-8) & (times <= C.N400_WIN[1] + 1e-8)
    cols = {k: np.full(len(tr), np.nan) for k in ("base_cp", "base_fz", "raw_cp", "raw_fz")}
    for section in C.SECTIONS:
        eeg = np.load(C.SHARED / subject / f"section_{section:03d}" / "eeg_data.npy", mmap_mode="r")
        if eeg.shape[0] != 65:
            raise SystemExit(f"{subject} section {section} eeg shape {eeg.shape}")
        # Pull only the four channels; mmap keeps the rest on disk.
        cp_sig = np.asarray(eeg[cp_ix], dtype=np.float64).mean(axis=0)
        fz_sig = np.asarray(eeg[fz_ix], dtype=np.float64)
        rows = np.where(tr["section"].to_numpy() == section)[0]
        for i in rows:
            start = onsets[(section, int(tr["word_index"].iat[i]))] - pre
            cp = cp_sig[start:start + n_t]
            fz = fz_sig[start:start + n_t]
            cols["base_cp"][i] = cp[base].mean()
            cols["base_fz"][i] = fz[base].mean()
            cols["raw_cp"][i] = cp[win].mean()
            cols["raw_fz"][i] = fz[win].mean()
    for k, v in cols.items():
        tr[k] = v
    # Stored dv came from float32 time courses; allow float32-level error.
    for roi in ("cp", "fz"):
        diff = (tr[f"raw_{roi}"] - tr[f"base_{roi}"]) - tr[f"dv_{roi}"]
        scale = np.maximum(1.0, tr[f"dv_{roi}"].abs())
        bad = ~(diff.abs() <= 1e-4 * scale)
        if bad.any():
            raise SystemExit(
                f"{subject} {roi}: {int(bad.sum())} rows do not reproduce dv "
                f"(max abs diff {float(diff.abs().max()):.3g})"
            )
    return tr[["subject", "section", "word_index", "item_id", "base_cp", "base_fz", "raw_cp", "raw_fz"]]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    words = pd.read_csv(C.STIM_DIR / "word_table.csv")
    onsets = {
        (int(r.section), int(r.word_index)): int(r.onset_sample)
        for r in words.itertuples(index=False)
    }
    montage = C.confirm_montage(65)
    cp_ix = [montage["index"][ch] for ch in C.CP_CHANNELS]
    fz_ix = montage["index"][C.FZ_CHANNEL]

    frames = []
    for s in C.load_subjects():
        f = subject_rows(s, onsets, cp_ix, fz_ix)
        frames.append(f)
        print(f"{s} rows {len(f)} reproduced dv", flush=True)
    bc = pd.concat(frames, ignore_index=True)

    lmm = pd.read_csv(C.STIM_DIR / "lmm_input.csv")
    n0 = len(lmm)
    lmm = lmm.merge(bc, on=["subject", "section", "word_index", "item_id"], how="left", validate="1:1")
    if len(lmm) != n0 or lmm[["base_cp", "raw_cp"]].isna().any().any():
        raise SystemExit("baseline merge lost or failed rows")
    lmm.to_csv(OUT / "lmm_input_basecov.csv", index=False)

    summary = {
        "n_rows": int(len(lmm)),
        "n_subjects": int(lmm["subject"].nunique()),
        "dv_reproduced": True,
        "baseline_s": list(C.BASELINE),
        "window_s": list(C.N400_WIN),
        "sd": {k: float(lmm[k].std()) for k in ("dv_cp", "raw_cp", "base_cp", "dv_fz", "raw_fz", "base_fz")},
        "corr_raw_base_cp": float(np.corrcoef(lmm["raw_cp"], lmm["base_cp"])[0, 1]),
        "corr_raw_base_fz": float(np.corrcoef(lmm["raw_fz"], lmm["base_fz"])[0, 1]),
    }
    (OUT / "extract_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
