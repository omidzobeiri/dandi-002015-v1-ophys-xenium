"""
build_multimodal_nwb.py
=======================

Make one multimodal NWB file (zarr backend, the same as the AIND processed
ophys NWB) for one ophys session of one Xenium mouse.

The file holds:
  * all ophys data of the session (copied from the processed NWB, all ROIs),
  * the stimulus (intervals; for the movie session also the movie frames in
    stimulus/templates),
  * running speed, facemap motion energy + motion SVD, LP-eye pupil / eye / CR,
  * cell types and cell x gene counts, but only for the ROIs that are
    coregistered to a Xenium cell in THIS session.

Coregistration comes from
scratch/ophys_transcriptomics_match_tables/<mouse>x_ophys_transcriptomics_match_all.csv.
In that table, 0 means "no match", and the ophys id is the ROI label, so the
NWB roi_table id = ophys id - 1. When a z-stack cell has more than one Xenium
match, the primary match is the lowest section number (the same rule as
build_long_match_tables_xenium.py). All matches go in a long table.

Stimulus times are corrected with the measured monitor delay from
scratch/ophys/natural_movies/timing_check/pd_delay.csv.

Usage (ophys_env; xenium_env does not have pynwb):
    python build_multimodal_nwb.py --mouse 786297 --session session_1
    python build_multimodal_nwb.py --mouse 786297 --session all
"""
import argparse
import glob
import json
import pickle
import re
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from hdmf.common import DynamicTable, ElementIdentifiers, VectorData
from hdmf_zarr import NWBZarrIO, ZarrDataIO
from numcodecs import Blosc
from pynwb import TimeSeries
from pynwb.image import ImageSeries
from pynwb.behavior import BehavioralTimeSeries, PupilTracking
from pynwb.epoch import TimeIntervals

REPO = Path("/root/capsule")
DATA = REPO / "data"
SCRATCH = REPO / "scratch"
MATCH_DIR = SCRATCH / "ophys_transcriptomics_match_tables"
TX_DIR = SCRATCH / "transcriptomics"
MOVIE_DIR = SCRATCH / "ophys" / "natural_movies" / "Natural_movie_files"
DELAY_CSV = SCRATCH / "ophys" / "natural_movies" / "timing_check" / "pd_delay.csv"
OUT_DIR = SCRATCH / "nwb_multimodal"

sys.path.insert(0, str(REPO / "code" / "Preprocessing"))
from analysis_utils import get_running_df, get_running_timestamps  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))

PLANES = [f"VISp_{i}" for i in range(8)]
CONTROL_PREFIXES = ("NegControl", "Unassigned", "Deprecated", "BLANK",
                    "Intergenic", "Genomic")
CELL_TYPE_LEVELS = ("class", "subclass", "supertype", "cluster")
FIRST_NON_NEURONAL_CLASS = 30  # ABC classes 30-34: Astro-Epen, OPC-Oligo, OEC, Vascular, Immune
COMPRESSOR = Blosc(cname="zstd", clevel=5, shuffle=Blosc.SHUFFLE)
DATA_INFO = pickle.load(open(REPO / "code/Preprocessing/Xenium/data_info_xenium.pkl", "rb"))


def zio(data, chunks=None):
    return ZarrDataIO(data=data, chunks=chunks, compressor=COMPRESSOR)


def newest(pattern):
    """Return the newest asset for a glob pattern. The asset name ends with its date."""
    hits = sorted(glob.glob(str(DATA / pattern)))
    return Path(hits[-1]) if hits else None


# ---------------------------------------------------------------- transcriptomics
def coregistered_rois(mouse, s_idx):
    """One row per ROI of session s_idx that has a z-stack match.

    Columns: plane, roi_id, czstack_id, unique_cell_id, n_xenium_matches,
    xenium_section, xenium_cell_id (primary; 0 if no Xenium match).
    Also returns the long table of all Xenium matches.
    """
    wide = pd.read_csv(MATCH_DIR / f"{mouse}x_ophys_transcriptomics_match_all.csv")
    wide = wide[wide.czstack_id > 0]
    xcols = sorted([c for c in wide if re.fullmatch(r"xenium_\d+_id", c)],
                   key=lambda c: int(c.split("_")[1]))
    rows, long = [], []
    for plane in PLANES:
        col = f"ophys_session_{s_idx}_{plane}_id"
        if col not in wide:
            continue
        sub = wide[wide[col] > 0]
        for _, r in sub.iterrows():
            hits = [(int(c.split("_")[1]), int(r[c])) for c in xcols if r[c] > 0]
            roi_id = int(r[col]) - 1
            rows.append(dict(plane=plane, roi_id=roi_id, czstack_id=int(r.czstack_id),
                             unique_cell_id=r.unique_cell_id, n_xenium_matches=len(hits),
                             xenium_section=hits[0][0] if hits else 0,
                             xenium_cell_id=hits[0][1] if hits else 0))
            for k, (sec, xid) in enumerate(hits):
                long.append(dict(plane=plane, roi_id=roi_id, unique_cell_id=r.unique_cell_id,
                                 xenium_section=sec, xenium_cell_id=xid, is_primary=k == 0))
    rois = pd.DataFrame(rows)
    assert not rois.duplicated(["plane", "roi_id"]).any(), "ROI listed twice in match table"
    return rois, pd.DataFrame(long)


def transcriptomics_tables(mouse, coreg):
    """Cell types and gene counts for the coregistered ROIs (primary match)."""
    key = ["xenium_section", "cell_id - xenium"]
    ct = pd.read_csv(TX_DIR / f"{mouse}x_cell_types.csv")
    gx = pd.read_csv(TX_DIR / f"{mouse}x_cellxgene.csv")
    # Genes = numeric count columns. Some tables also have text id columns (for example 'index').
    genes = [c for c in gx.columns
             if c not in ("mouse_id", *key) and not c.startswith(CONTROL_PREFIXES)
             and pd.api.types.is_numeric_dtype(gx[c])]
    left = coreg.rename(columns={"xenium_cell_id": "cell_id - xenium"})
    ct_m = left.merge(ct.drop(columns="mouse_id"), on=key, how="left")
    gx_m = left.merge(gx[key + genes], on=key, how="left")
    assert len(ct_m) == len(coreg) and len(gx_m) == len(coreg)
    return ct_m, gx_m[genes].to_numpy(), genes


def add_transcriptomics(nwb, mouse, s_idx, rois, long):
    coreg = rois[rois.xenium_cell_id > 0].sort_values(["plane", "roi_id"]).reset_index(drop=True)
    ct, counts, genes = transcriptomics_tables(mouse, coreg)
    # Remove non-neuronal cells (ABC classes 30-34). Cells with no cell-type mapping stay.
    class_num = ct["class_name"].str.extract(r"^(\d+)")[0].astype(float)
    keep = ~(class_num >= FIRST_NON_NEURONAL_CLASS).to_numpy()
    n_non_neuronal = int((~keep).sum())
    coreg = coreg[keep].reset_index(drop=True)
    ct = ct[keep].reset_index(drop=True)
    counts = counts[keep]
    n_missing_ct = int(ct["class_name"].isna().sum())
    n_missing_gx = int(np.isnan(counts).any(axis=1).sum())
    counts = np.nan_to_num(counts, nan=-1).astype(np.int32)

    xen = DATA_INFO[int(mouse)]["xenium"]
    mod = nwb.create_processing_module(
        "transcriptomics",
        f"Xenium spatial transcriptomics for ROIs coregistered in this session "
        f"(ophys_session_{s_idx} of the match table). Sources: {xen['processed']}, "
        f"{xen['cell_types']}, {MATCH_DIR.name}/{mouse}x_ophys_transcriptomics_match_all.csv. "
        f"Primary Xenium match = lowest section number. Only neurons: {n_non_neuronal} ROIs "
        f"matched to non-neuronal cells (ABC classes 30-34) are not included. Rows with no "
        f"cell type: {n_missing_ct}; rows with no gene counts (counts = -1): {n_missing_gx}.")

    id_cols = [
        ("plane", "Imaging plane (processing module name)", coreg.plane.to_numpy()),
        ("roi_id", "Row id in processing/<plane>/image_segmentation/roi_table", coreg.roi_id.to_numpy()),
        ("unique_cell_id", "Cell id that is the same in all sessions of this mouse", coreg.unique_cell_id.to_numpy()),
        ("czstack_id", "Cell id in the cortical z-stack", coreg.czstack_id.to_numpy()),
        ("xenium_section", "Xenium section of the primary match", coreg.xenium_section.to_numpy()),
        ("xenium_cell_id", "Xenium cell id of the primary match", coreg.xenium_cell_id.to_numpy()),
        ("n_xenium_matches", "Number of Xenium cells matched to this z-stack cell", coreg.n_xenium_matches.to_numpy()),
    ]

    def table(name, desc, extra):
        cols = [VectorData(name=n, description=d, data=v) for n, d, v in id_cols + extra]
        return DynamicTable(name=name, description=desc, columns=cols)

    ct_cols = []
    for lvl in CELL_TYPE_LEVELS:
        ct_cols += [(f"{lvl}_label", f"ABC atlas {lvl} label", ct[f"{lvl}_label"].fillna("").to_numpy()),
                    (f"{lvl}_name", f"ABC atlas {lvl} name", ct[f"{lvl}_name"].fillna("").to_numpy()),
                    (f"{lvl}_bootstrapping_probability", f"Mapping bootstrapping probability, {lvl}",
                     ct[f"{lvl}_bootstrapping_probability"].to_numpy(float))]
    ct_cols.append(("cluster_alias", "ABC atlas cluster alias",
                    ct["cluster_alias"].fillna(-1).astype(int).to_numpy()))
    mod.add(table("cell_types", "Cell-type mapping (ABC atlas taxonomy CS20230722) of coregistered ROIs", ct_cols))

    mod.add(table("cell_by_gene", "Xenium transcript counts of coregistered ROIs. Column order of "
                  "'counts' = rows of the 'genes' table. Control codewords are not included.",
                  [("counts", "Transcript counts, shape (n_cells, n_genes)", zio(counts))]))

    mod.add(DynamicTable(name="genes", description="Genes in the Xenium panel, in the order of cell_by_gene.counts",
                         columns=[VectorData(name="gene_name", description="Gene or probe name",
                                             data=np.array(genes))]))

    lg = long.merge(coreg[["plane", "roi_id"]], on=["plane", "roi_id"]) if len(long) else long
    mod.add(DynamicTable(
        name="xenium_matches", description="All Xenium matches of the coregistered ROIs (one row per match)",
        columns=[VectorData(name=c, description=d, data=lg[c].to_numpy()) for c, d in [
            ("plane", "Imaging plane"), ("roi_id", "Row id in the plane roi_table"),
            ("unique_cell_id", "Cross-session cell id"), ("xenium_section", "Xenium section"),
            ("xenium_cell_id", "Xenium cell id"), ("is_primary", "True for the primary match")]]))
    return coreg, n_missing_ct, n_missing_gx, n_non_neuronal


# ---------------------------------------------------------------- ROI table columns
def add_roi_columns(nwb, rois, coreg):
    for plane in PLANES:
        if plane not in nwb.processing:
            continue
        rt = nwb.processing[plane]["image_segmentation"]["roi_table"]
        ids = np.asarray(rt.id[:])
        r = rois[rois.plane == plane].set_index("roi_id")
        c = coreg[coreg.plane == plane].reset_index().set_index("roi_id")["index"]
        rt.add_column("unique_cell_id", "Cross-session cell id ('' if no z-stack match)",
                      data=[str(r.unique_cell_id.get(i, "")) for i in ids])
        rt.add_column("czstack_id", "Cell id in the cortical z-stack (0 if no match)",
                      data=np.array([int(r.czstack_id.get(i, 0)) for i in ids]))
        rt.add_column("is_coregistered", "True if the ROI has a row in processing/transcriptomics: "
                      "a Xenium match in this session that is not a non-neuronal cell",
                      data=np.array([i in c.index for i in ids]))
        rt.add_column("transcriptomics_row", "Row in processing/transcriptomics tables (-1 if none)",
                      data=np.array([int(c.get(i, -1)) for i in ids]))


# ---------------------------------------------------------------- stimulus
def delay_shift_s(mouse, raw):
    d = pd.read_csv(DELAY_CSV)
    date = re.search(r"\d{4}-\d{2}-\d{2}", raw).group(0)
    row = d[(d["mouse"].astype(str) == str(mouse)) & (d["date"] == date)]
    assert len(row) == 1, f"pd_delay.csv: {len(row)} rows for {mouse} {date}"
    row = row.iloc[0]
    return (row.measured_delay_ms - row.code_delay_ms) / 1000.0, row


def add_stimulus(nwb, mouse, raw, is_movie):
    st = pd.read_csv(glob.glob(str(DATA / raw / "behavior" / "*_stim_table.csv"))[0])
    shift, row = delay_shift_s(mouse, raw)
    note = (f"Times = stim-table time - {row.code_delay_ms:.1f} ms + measured monitor delay "
            f"{row.measured_delay_ms:.1f} ms (shift {shift*1000:+.1f} ms). Original times are in "
            f"*_uncorrected columns.")

    def intervals(name, desc, df, cols):
        cols = [("start_time", "Start time (s), corrected", df.start_time.to_numpy() + shift),
                ("stop_time", "Stop time (s), corrected", df.stop_time.to_numpy() + shift)] + \
               [(c, d, df[c].to_numpy()) for c, d in cols] + \
               [("start_time_uncorrected", "Start time in the stim table (s)", df.start_time.to_numpy()),
                ("stop_time_uncorrected", "Stop time in the stim table (s)", df.stop_time.to_numpy())]
        nwb.add_time_intervals(TimeIntervals(
            name=name, description=f"{desc} {note}",
            id=ElementIdentifiers(name="id", data=np.arange(len(df))),
            columns=[VectorData(name=c, description=d, data=v) for c, d, v in cols]))

    sp = st[st.stim_name == "spontaneous"]
    intervals("spontaneous_presentations", "Gray-screen spontaneous blocks.", sp, [])
    ev = st[st.stim_name != "spontaneous"].copy()
    ev["stim_block"] = ev.stim_block.astype(int)
    ev["stim_index"] = ev.stim_index.astype(int)
    if not is_movie:
        cols = [("stim_name", "Stimulus name"), ("stim_block", "Stimulus block"),
                ("contrast", "Grating contrast"), ("temporal_frequency", "Temporal frequency (Hz)"),
                ("spatial_frequency", "Spatial frequency (cycles/deg)"),
                ("orientation", "Orientation (deg)"), ("stim_index", "Stimulus index")]
        intervals("drifting_gratings_presentations", "Drifting-grating trials.", ev, cols)
        return len(ev), shift

    ev["frame"] = ev.frame.astype(int)
    ev = ev.rename(columns={"stim_name": "movie_name"})
    intervals("natural_movie_presentations", "One row per displayed movie frame.", ev,
              [("movie_name", "Movie name = name of the ImageSeries in stimulus/templates"),
               ("frame", "Frame index into the movie"), ("stim_block", "Stimulus block"),
               ("stim_index", "Stimulus index")])

    for name in sorted(ev.movie_name.unique()):
        arr = np.load(MOVIE_DIR / f"ds_warped_15_{name}.npy", mmap_mode="r")
        n_used = int(ev.loc[ev.movie_name == name, "frame"].max()) + 1
        assert n_used <= arr.shape[0], f"{name}: stim table uses {n_used} frames, file has {arr.shape[0]}"
        nwb.add_stimulus_template(ImageSeries(
            name=name, data=zio(arr, chunks=(30, arr.shape[1], arr.shape[2])), unit="n.a.",
            rate=30.0, starting_time=0.0, format="raw",
            description=f"Warped natural movie frames (uint8, frames x height x width) from "
                        f"{MOVIE_DIR.relative_to(REPO)}/ds_warped_15_{name}.npy. Frame i is shown "
                        f"in natural_movie_presentations rows with movie_name='{name}', frame=i."))
    return len(ev), shift


# ---------------------------------------------------------------- behavior
def clean_ts(df):
    df = df[np.isfinite(df.timestamps)]
    assert np.all(np.diff(df.timestamps.to_numpy()) > 0), "timestamps not increasing"
    return df


def add_behavior(nwb, raw):
    fm = newest(f"{raw}_facemap-motion-svd_*")
    lp = newest(f"{raw}_lp-eye_*")
    fm_txt = fm.name if fm else "no facemap asset for this session"
    lp_txt = lp.name if lp else "no LP-eye asset for this session"
    mod = nwb.create_processing_module(
        "behavior", f"Running (sync + stim.pkl), facemap ({fm_txt}), LP-eye ({lp_txt}). "
                    f"Timestamps on the sync clock (s).")
    info = dict(facemap_asset=fm.name if fm else None, lp_eye_asset=lp.name if lp else None)

    beh = DATA / raw / "behavior"
    ts = get_running_timestamps(glob.glob(str(beh / "*.h5"))[0])
    run = get_running_df(glob.glob(str(beh / "*.pkl"))[0], ts)
    run = clean_ts(run)
    rs = TimeSeries(name="running_speed", data=zio(run.speed.to_numpy(np.float32)), unit="cm/s",
                    timestamps=zio(run.timestamps.to_numpy()),
                    description="Wheel running speed, low-pass filtered (analysis_utils.get_running_df)")
    rsf = TimeSeries(name="running_speed_filtered", data=zio(run.speed_filtered.to_numpy(np.float32)),
                     unit="cm/s", timestamps=rs,
                     description="Rolling mean (50 samples) of |running_speed|")
    mod.add(BehavioralTimeSeries(name="running", time_series=[rs, rsf]))
    info["running_samples"] = len(run)

    if fm is not None:
        motion, svd = [], []
        with pd.HDFStore(fm / "facemap_table.h5", "r") as store:
            regions = [r for r in ("face", "nose", "behavior") if f"/{r}/motion" in store.keys()]
        info["facemap_regions"] = regions
        for region in regions:
            m = clean_ts(pd.read_hdf(fm / "facemap_table.h5", f"/{region}/motion"))
            v = clean_ts(pd.read_hdf(fm / "facemap_table.h5", f"/{region}/motsvd"))
            assert np.array_equal(m.timestamps.to_numpy(), v.timestamps.to_numpy())
            t0 = TimeSeries(name=f"{region}_motion_energy", data=zio(m.motion_energy.to_numpy(np.float32)),
                            unit="a.u.", timestamps=zio(m.timestamps.to_numpy()),
                            description=f"Facemap motion energy, {region} ROI")
            motion += [t0,
                       TimeSeries(name=f"{region}_motion_energy_clean", unit="a.u.", timestamps=t0,
                                  data=zio(m.motion_energy_clean.to_numpy(np.float32)),
                                  description=f"Motion energy with keyframe artifacts removed, {region} ROI"),
                       TimeSeries(name=f"{region}_is_keyframe_contaminated", unit="n.a.", timestamps=t0,
                                  data=zio(m.is_keyframe_contaminated.to_numpy().astype(np.uint8)),
                                  description="1 = frame contaminated by a video keyframe")]
            pcs = [c for c in v.columns if c.startswith("motsvd_")]
            svd.append(TimeSeries(name=f"{region}_motion_svd", unit="a.u.", timestamps=t0,
                                  data=zio(v[pcs].to_numpy(np.float32), chunks=(10000, len(pcs))),
                                  description=f"Facemap motion SVD, {region} ROI, shape (frames, {len(pcs)})"))
        mod.add(BehavioralTimeSeries(name="facemap_motion_energy", time_series=motion))
        mod.add(BehavioralTimeSeries(name="facemap_motion_svd", time_series=svd))

    if lp is None:
        return info
    e = clean_ts(pd.read_hdf(lp / "eye_tracking_table.h5", "/eye_tracking"))
    t_eye = None
    groups = {}
    for part in ("pupil", "eye", "cr"):
        series = []
        for q in ("area", "width", "height", "center_x", "center_y", "phi",
                  "average_confidence", "is_bad_frame"):
            col = f"{part}_{q}"
            data = e[col].to_numpy(np.uint8 if q == "is_bad_frame" else np.float32)
            unit = {"area": "pixels^2", "phi": "radians", "average_confidence": "n.a.",
                    "is_bad_frame": "n.a."}.get(q, "pixels")
            kw = dict(timestamps=t_eye) if t_eye is not None else dict(timestamps=zio(e.timestamps.to_numpy()))
            ts_ = TimeSeries(name=col, data=zio(data), unit=unit, description=f"LP-eye ellipse fit: {col}", **kw)
            if t_eye is None:
                t_eye = ts_
            series.append(ts_)
        groups[part] = series
    mod.add(PupilTracking(name="pupil_tracking", time_series=groups["pupil"]))
    mod.add(BehavioralTimeSeries(name="eye_tracking", time_series=groups["eye"]))
    mod.add(BehavioralTimeSeries(name="corneal_reflection_tracking", time_series=groups["cr"]))
    info["pupil_bad_frac"] = float(e.pupil_is_bad_frame.mean())
    return info


# ---------------------------------------------------------------- ophys source
def aind_metadata_ok(path, mouse, raw):
    """True if the AIND NWB has the real subject, session id and plane metadata.

    Many AIND NWB files have placeholder metadata (subject 'subject', session_id
    'ophys_session', start time = processing time, empty lab metadata)."""
    with NWBZarrIO(str(path), "r") as io:
        n = io.read()
        lm = list(n.lab_meta_data.values())
        return (n.subject is not None and n.subject.subject_id == str(mouse)
                and n.session_id == raw and bool(lm) and bool(lm[0].fields))


def ophys_source(mouse, raw, proc):
    """Return (path, kind) of the ophys NWB to copy from.

    kind = 'aind' (the AIND NWB, when its metadata is correct) or 'converted'
    (made by make_ophys_nwb.py from the per-plane files; the data are identical)."""
    hits = glob.glob(str(DATA / proc / "*.nwb"))
    if len(hits) == 1 and aind_metadata_ok(hits[0], mouse, raw):
        return Path(hits[0]), "aind"
    base = OUT_DIR / "ophys_base" / f"{proc}.nwb"
    done = base.with_suffix(".done")
    if not done.exists():
        from make_ophys_nwb import make
        shutil.rmtree(base, ignore_errors=True)
        base.parent.mkdir(parents=True, exist_ok=True)
        print(f"  making ophys base NWB {base.name} from per-plane files", flush=True)
        make(DATA / proc, DATA / raw, base)
        done.write_text("ok")
    return base, "converted"


# ---------------------------------------------------------------- main
def build(mouse, session):
    s_idx = int(session.split("_")[1])
    sinfo = DATA_INFO[int(mouse)]["ophys"][session]
    raw, proc = sinfo["raw"], sinfo["processed"]
    is_movie = sinfo["session_type"] == "STAGE_0"
    src, src_kind = ophys_source(mouse, raw, proc)
    out_dir = OUT_DIR / str(mouse)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{raw}_multimodal.nwb"
    assert not out.exists(), f"{out} exists; remove it first"
    t0 = time.time()
    print(f"[{mouse} {session}] {raw} ({'movie' if is_movie else 'DG'})", flush=True)

    rois, long = coregistered_rois(mouse, s_idx)
    with NWBZarrIO(str(src), "r", load_namespaces=True) as io:
        nwb = io.read()
        coreg, n_mct, n_mgx, n_nn = add_transcriptomics(nwb, mouse, s_idx, rois, long)
        add_roi_columns(nwb, rois, coreg)
        n_stim, shift = add_stimulus(nwb, mouse, raw, is_movie)
        binfo = add_behavior(nwb, raw)
        print(f"  built in memory ({time.time()-t0:.0f} s); writing {out}", flush=True)
        with NWBZarrIO(str(out), "w") as eio:
            eio.export(src_io=io, nwbfile=nwb, write_args=dict(link_data=False))

    summary = dict(mouse=str(mouse), session=session, raw=raw, processed=proc, output=str(out),
                   ophys_source=src_kind, is_movie=is_movie, n_rois_zstack_matched=int(len(rois)),
                   n_coregistered=int(len(coreg)), n_non_neuronal_removed=n_nn,
                   n_multi_xenium=int((coreg.n_xenium_matches > 1).sum()),
                   n_missing_cell_type=n_mct, n_missing_counts=n_mgx,
                   n_stim_rows=int(n_stim), stim_shift_ms=round(shift * 1000, 2),
                   seconds=round(time.time() - t0), **binfo)
    json.dump(summary, open(out_dir / f"{raw}_multimodal_summary.json", "w"), indent=2)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mouse", required=True)
    ap.add_argument("--session", required=True, help="session_0..session_3 or all")
    a = ap.parse_args()
    sessions = [f"session_{i}" for i in range(4)] if a.session == "all" else [a.session]
    for s in sessions:
        build(a.mouse, s)
