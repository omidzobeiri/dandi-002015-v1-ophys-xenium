"""
Read the parts of one DANDI:002015 NWB file (zarr backend).

    from v1_ophys_xenium import open_nwb, roi_table, traces, transcriptomics

    with open_nwb(path) as nwb:
        rois = roi_table(nwb, "VISp_1")
        t, dff, roi_ids = traces(nwb, "VISp_1", "dff")
        tx = transcriptomics(nwb)

The trace data are read from disk only when you index them (for example
dff[:, 0]), so keep the file open while you read data.
"""
from contextlib import contextmanager

import numpy as np
import pandas as pd
from hdmf_zarr import NWBZarrIO

# kind -> (data interface, series name or None if the interface is the series)
TRACE_KINDS = {
    "dff": ("dff_timeseries", "dff_timeseries"),
    "events": ("event_timeseries", None),
    "raw": ("raw_timeseries", "ROI_fluorescence_timeseries"),
    "neuropil": ("neuropil_fluorescence_timeseries", None),
    "corrected": ("neuropil_corrected_timeseries", None),
}


@contextmanager
def open_nwb(path):
    """Open an NWB zarr file for reading. Use it in a 'with' statement."""
    io = NWBZarrIO(str(path), "r")
    try:
        yield io.read()
    finally:
        io.close()


# ---------------------------------------------------------------- ophys
def imaging_planes(nwb):
    """Names of the imaging-plane processing modules (VISp_0 ... VISp_7)."""
    return sorted(p for p in nwb.processing if p.startswith("VISp_"))


def plane_depths(nwb):
    """Imaging depth (µm) of each plane. The plane number is not in depth order."""
    return {p: int(nwb.lab_meta_data[p].imaging_depth) for p in imaging_planes(nwb)
            if p in nwb.lab_meta_data}


def roi_table(nwb, plane):
    """ROI table of one plane as a DataFrame, without the large image_mask column.

    The builder added: unique_cell_id, czstack_id, is_coregistered, transcriptomics_row.
    """
    rt = nwb.processing[plane]["image_segmentation"]["roi_table"]
    cols = [c for c in rt.colnames if c != "image_mask"]
    return rt.to_dataframe(exclude={"image_mask"})[cols]


def plane_summary(nwb):
    """One row per plane: depth, number of ROIs, soma ROIs, z-stack matches and coregistered ROIs."""
    depth = plane_depths(nwb)
    rows = []
    for p in imaging_planes(nwb):
        t = roi_table(nwb, p)
        rows.append(dict(plane=p, depth_um=depth.get(p), n_rois=len(t), n_soma=int(t.is_soma.sum()),
                         n_zstack_matched=int((t.czstack_id > 0).sum()),
                         n_coregistered=int(t.is_coregistered.sum())))
    return pd.DataFrame(rows)


def traces(nwb, plane, kind="dff"):
    """Return (timestamps, data [time x roi], roi ids) for one plane.

    kind: 'dff', 'events', 'raw', 'neuropil' or 'corrected'.
    data is a lazy zarr array: index it to read values.
    """
    iface, name = TRACE_KINDS[kind]
    ts = nwb.processing[plane][iface]
    if name is not None:
        ts = ts[name]
    return ts.timestamps[:], ts.data, np.asarray(ts.rois.data[:])


# ---------------------------------------------------------------- transcriptomics
def transcriptomics(nwb):
    """Return a dict of DataFrames from processing/transcriptomics.

    cell_types:     one row per coregistered ROI (plane, roi_id, ids, class ... cluster labels)
    counts:         cell x gene counts, same rows as cell_types, one column per gene
    xenium_matches: all Xenium matches; is_primary marks the match used in the other tables
    """
    tx = nwb.processing["transcriptomics"]
    cell_types = tx["cell_types"].to_dataframe()
    genes = np.asarray(tx["genes"]["gene_name"].data[:])
    counts = pd.DataFrame(tx["cell_by_gene"]["counts"].data[:], columns=genes, index=cell_types.index)
    return dict(cell_types=cell_types, counts=counts,
                xenium_matches=tx["xenium_matches"].to_dataframe())


def celltype_traces(nwb, cell_types, subclass, kind="dff"):
    """Return (timestamps, [time x cells] array, cell_types rows) for all coregistered ROIs of a subclass.

    The planes have slightly different frame times, so all planes are put on the
    clock of the first plane with linear interpolation.
    """
    rows = cell_types[cell_types.subclass_name == subclass]
    out, t_ref = [], None
    for plane, grp in rows.groupby("plane"):
        t_p, data, ids = traces(nwb, plane, kind)
        col = {r: i for i, r in enumerate(ids)}
        x = data[:, [col[r] for r in grp.roi_id]]
        if t_ref is None:
            t_ref = t_p
        out.append(np.column_stack([np.interp(t_ref, t_p, x[:, j]) for j in range(x.shape[1])]))
    return t_ref, np.hstack(out), rows


def coregistered_cells(path):
    """The coregistered ROIs of one file: plane, roi_id, unique_cell_id, subclass_name."""
    with open_nwb(path) as nwb:
        return nwb.processing["transcriptomics"]["cell_types"].to_dataframe()[
            ["plane", "roi_id", "unique_cell_id", "subclass_name"]]


# ---------------------------------------------------------------- stimulus
def is_movie_session(nwb):
    return "natural_movie_presentations" in nwb.intervals


def stimulus_tables(nwb):
    """Return (stimulus table, spontaneous table) as DataFrames.

    The stimulus table is natural_movie_presentations (one row per shown movie
    frame) or drifting_gratings_presentations (one row per trial).
    """
    name = "natural_movie_presentations" if is_movie_session(nwb) else "drifting_gratings_presentations"
    return nwb.intervals[name].to_dataframe(), nwb.intervals["spontaneous_presentations"].to_dataframe()


# ---------------------------------------------------------------- behavior
def behavior_series(nwb, interface, name):
    """Return (timestamps, lazy data) of one series in processing/behavior.

    Example: behavior_series(nwb, "pupil_tracking", "pupil_area").
    """
    s = nwb.processing["behavior"][interface][name]
    return s.timestamps[:], s.data


def align_to(t_target, t_source, x):
    """Linear interpolation of a 1-D signal x(t_source) onto t_target. NaN samples are ignored."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    return np.interp(t_target, t_source[ok], x[ok])
