# DANDI:002015 — mouse V1 two-photon imaging with Xenium spatial transcriptomics of the same neurons

This repository goes with the dataset [DANDI:002015](https://dandiarchive.org/dandiset/002015).
It gives:

- `v1_ophys_xenium/`: a small Python package to find, download and read the NWB files.
- `notebooks/load_multimodal_nwb.ipynb`: a notebook that shows each part of one NWB file.
- `build_scripts/`: the scripts that made the NWB files (for reference only; see [Build scripts](#build-scripts)).

> **Access.** The dandiset is embargoed until **2027-12-28**. Until that date, only users that have access to the dandiset can download the files. See [Get access and an API key](#2-get-access-and-an-api-key).

## Contents

1. [The dataset](#the-dataset)
2. [The NWB file layout](#the-nwb-file-layout)
3. [Get started](#get-started)
4. [Use the helper functions](#use-the-helper-functions)
5. [Use the notebook](#use-the-notebook)
6. [Important notes](#important-notes)
7. [Build scripts](#build-scripts)
8. [Citation and contact](#citation-and-contact)

## The dataset

The dataset connects the activity of neurons in mouse primary visual cortex (VISp) to their molecular cell type.

| Item | Value |
|---|---|
| Mice | 9 (Snap25-IRES2-Cre; Oi4 jGCaMP8s, pan-neuronal) |
| Sessions | 4 per mouse, 36 in total, one NWB file per session |
| Size | 127.8 GB (about 9 GB for a movie session, 1.4–2.2 GB for a grating session) |
| Imaging | Multiplane two-photon mesoscope, 8 planes in VISp, about 10.7 Hz per plane |
| Transcriptomics | 10x Xenium, 299-gene panel, cells mapped to the Allen Brain Cell (ABC) Atlas taxonomy |
| Coregistered neurons | 959–2,127 ROIs per session (57,482 ROI-session pairs in total) |
| Format | NWB 2, zarr backend (`.nwb.zarr`), the same format as the AIND processed ophys NWB |

**Experiment.** We recorded the activity of neurons in awake, head-fixed mice. Each mouse had 4 imaging sessions:

- 1 natural-movie session: 55 movie clips and 1 repeated test clip (56 movies).
- 3 drifting-grating sessions: 2 s gratings in a contrast series and in a temporal-frequency series (1,120 trials per session).
- All sessions also have spontaneous (gray-screen) blocks.

During each session we also recorded:

- wheel running speed,
- pupil, eye and corneal-reflection position and size (Lightning Pose keypoints with ellipse fits),
- face, nose and body motion (Facemap motion energy and motion SVD).

**Transcriptomics.** After the imaging, we measured gene expression in serial tissue sections with 10x Xenium. We registered the Xenium cells to an in vivo cortical z-stack, and the z-stack to the imaged ROIs. We mapped the Xenium cells to the ABC Atlas taxonomy (class, subclass, supertype and cluster).

**Files.** The file names follow the DANDI rules:

```
sub-<mouse>/sub-<mouse>_ses-<YYYYMMDDTHHMMSS>_behavior+image+ophys.nwb.zarr   natural-movie session (has the movie frames)
sub-<mouse>/sub-<mouse>_ses-<YYYYMMDDTHHMMSS>_behavior+ophys.nwb.zarr         drifting-grating session
```

Mouse ids: 778174, 786297, 797371, 810976, 816460, 816965, 818629, 827542, 832701.

## The NWB file layout

Each file holds one imaging session:

```
/processing/VISp_0 ... VISp_7          one module per imaging plane (all ROIs)
    dff_timeseries/dff_timeseries      dF/F            (time x ROI)
    event_timeseries                   detected events (time x ROI)
    raw_timeseries/ROI_fluorescence_timeseries
    neuropil_fluorescence_timeseries
    neuropil_corrected_timeseries
    images                             average / max projection, segmentation_mask_image
    image_segmentation/roi_table       ROI masks and ROI metadata (+ 4 cross-modal columns)
/processing/transcriptomics            only for ROIs coregistered in THIS session, only neurons
    cell_types                         one row per coregistered ROI: ids + class/subclass/supertype/cluster
    cell_by_gene                       counts (n_cells x n_genes), same rows as cell_types
    genes                              gene names, in the column order of counts
    xenium_matches                     all Xenium matches; is_primary marks the one used above
/processing/behavior
    running                            running_speed (cm/s), running_speed_filtered
    pupil_tracking / eye_tracking / corneal_reflection_tracking
    facemap_motion_energy              <region>_motion_energy(_clean), _is_keyframe_contaminated
    facemap_motion_svd                 <region>_motion_svd (frames x 100 PCs)
/intervals
    drifting_gratings_presentations    grating sessions: one row per trial
    natural_movie_presentations        movie session: one row per shown movie frame
    spontaneous_presentations          gray-screen blocks
/stimulus/templates                    movie session: one ImageSeries per movie (uint8, frames x 300 x 480)
/general/lab_meta_data/VISp_<k>        imaging depth of each plane
```

Each `roi_table` has 4 cross-modal columns:

| Column | Meaning |
|---|---|
| `unique_cell_id` | Cell id that is the same in the 4 sessions of the mouse (`''` if the ROI has no z-stack match) |
| `czstack_id` | Cell id in the cortical z-stack (0 if no match) |
| `is_coregistered` | True if the ROI has a row in the transcriptomics tables of this file |
| `transcriptomics_row` | Row in the `processing/transcriptomics` tables (−1 if none) |

All times (ophys, behavior and stimulus) are in seconds on one sync clock.

## Get started

### 1. Make the environment

You need Python 3.10 or later. Clone the repository, then make the environment with conda **or** with pip.

With conda:

```bash
git clone https://github.com/omidzobeiri/dandi-002015-v1-ophys-xenium.git
cd dandi-002015-v1-ophys-xenium
conda env create -f environment.yml
conda activate dandi002015
```

With pip:

```bash
git clone https://github.com/omidzobeiri/dandi-002015-v1-ophys-xenium.git
cd dandi-002015-v1-ophys-xenium
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

The main packages are `pynwb`, `hdmf-zarr` and `zarr` 2 (to read the files), `dandi` (to download the files), and `pandas`, `numpy` and `matplotlib`. We tested the package and the notebook with Python 3.12, pynwb 4.2, hdmf 6.2, hdmf-zarr 0.13, zarr 2.18 and dandi 0.80.

### 2. Get access and an API key

While the dandiset is embargoed:

1. Make an account on [dandiarchive.org](https://dandiarchive.org).
2. Ask the contact person (see [Citation and contact](#citation-and-contact)) to add your account to the dandiset.
3. Copy your API key: log in, click your initials at the top right, and copy the key.
4. Set the key in your shell before you download:

```bash
export DANDI_API_KEY=<your key>
```

After the embargo ends, the files are public and you do not need a key.

### 3. Download the files

The files are large, so download only the sessions that you need. With Python:

```python
import v1_ophys_xenium as vx

sessions = vx.list_sessions()                     # all 36 files: subject, session, session_type, path, size_gb
print(sessions)

vx.download_sessions("sub-786297", "data")        # all 4 sessions of one mouse (about 16 GB)
vx.download_sessions(sessions.path[5], "data")    # or one file
```

Or with the `dandi` command:

```bash
dandi download "dandi://dandi/002015@draft/sub-786297/" -o data
```

A download that stops can be started again. The files that are complete are skipped.

## Use the helper functions

The package `v1_ophys_xenium` has these functions:

| Function | What it does |
|---|---|
| `list_sessions(subject=None)` | Table of the files in the dandiset |
| `download_sessions(paths, out_dir)` | Download files or a whole `sub-<mouse>` folder |
| `local_sessions(root)` | Table of the files that are on your disk |
| `open_nwb(path)` | Open a file (use it in a `with` statement) |
| `imaging_planes(nwb)`, `plane_depths(nwb)` | Plane names and imaging depths (µm) |
| `plane_summary(nwb)` | ROIs, soma ROIs, z-stack matches and coregistered ROIs per plane |
| `roi_table(nwb, plane)` | ROI table of one plane (without the large image masks) |
| `traces(nwb, plane, kind)` | Timestamps, lazy data (time × ROI) and ROI ids. `kind` is `dff`, `events`, `raw`, `neuropil` or `corrected` |
| `transcriptomics(nwb)` | `cell_types`, `counts` (cell × gene) and `xenium_matches` as DataFrames |
| `celltype_traces(nwb, cell_types, subclass)` | Traces of all coregistered ROIs of one subclass, on one clock |
| `stimulus_tables(nwb)` | The stimulus table and the spontaneous table |
| `is_movie_session(nwb)` | True for the natural-movie session |
| `behavior_series(nwb, interface, name)` | Timestamps and lazy data of one behavior series |
| `align_to(t_target, t_source, x)` | Put a signal on other timestamps (for example behavior on ophys frames) |
| `coregistered_cells(path)` | Coregistered ROIs of one file, to join sessions on `unique_cell_id` |

Example: the dF/F of all Pvalb neurons in one session, and the running speed on the same time base.

```python
import v1_ophys_xenium as vx

files = vx.local_sessions("data")
path = files[(files.subject == "786297") & (files.session_type == "drifting_gratings")].path.iloc[0]

with vx.open_nwb(path) as nwb:
    tx = vx.transcriptomics(nwb)
    t, pvalb, rows = vx.celltype_traces(nwb, tx["cell_types"], "052 Pvalb Gaba")
    t_run, run = vx.behavior_series(nwb, "running", "running_speed")
    run_on_ophys = vx.align_to(t, t_run, run[:])
    grating_trials, spont = vx.stimulus_tables(nwb)

print(pvalb.shape, len(rows), "Pvalb neurons")
```

Example: the same neurons in two sessions.

```python
a = vx.coregistered_cells(files.path.iloc[0])
b = vx.coregistered_cells(files.path.iloc[1])
both = a.merge(b, on="unique_cell_id", suffixes=("_a", "_b"))
```

## Use the notebook

```bash
conda activate dandi002015          # or: source .venv/bin/activate
jupyter lab notebooks/load_multimodal_nwb.ipynb
```

1. In the first code cell, set `DATA_DIR` (default `../data`, the `data` folder in the repository), `MOUSE` and `SESSION`.
2. If the files of `MOUSE` are not in `DATA_DIR`, the cell downloads them (set `DANDI_API_KEY` first).
3. Run all cells.

The notebook has 7 parts: open the file; ophys data (Figure 1: field of view and coregistered ROIs); transcriptomics (Figure 2: subclasses, Figure 3: marker genes); stimulus (Figure 4: responses by subclass); behavior (Figure 5); natural movies (Figure 6); and the same cell in other sessions.
The figure captions give the results for the default file (mouse 786297, session 2025-05-13).

## Important notes

- **Coregistered set per session.** The transcriptomics tables of a file hold only the ROIs that are coregistered in that session. The set is different in each session. To find a neuron in other sessions, join on `unique_cell_id`.
- **Only neurons.** ROIs matched to non-neuronal Xenium cells (ABC classes 30–34) are not in the tables (2,869 ROI matches in all files). The description of `processing/transcriptomics` gives the number for each file.
- **Unmapped ROIs.** Some ROIs have a Xenium match but no cell-type mapping. They stay in the tables, with names `''` and probabilities NaN.
- **ROI ids.** `cell_types.roi_id` is the id in the plane's `roi_table`. In `segmentation_mask_image`, label *k* is `roi_table` row *k* − 1.
- **Plane order.** The plane number (`VISp_0` … `VISp_7`) is not in depth order. Use `plane_depths(nwb)`.
- **Stimulus times.** The stimulus times are corrected with the monitor delay that we measured with a photodiode in each session (a shift from −0.9 to +18.6 ms, different for each session). The original times are in the `*_uncorrected` columns. The table description gives the shift.
- **Missing behavior data.** Mice 778174 and 797371 have no pupil and no Facemap data (running speed only). The movie session of mouse 810976 has no nose region in Facemap.
- **Bad pupil frames.** Use `pupil_is_bad_frame` (and the same series for eye and CR) to remove bad frames.
- **Lazy data.** The trace data are read only when you index them. Keep the file open while you read.

## Build scripts

The folder `build_scripts/` holds the scripts that made the files. They use internal paths and source data that are not in the dandiset, so you cannot run them as they are. We include them to show how each value in the files was made.

- `build_multimodal_nwb.py`: makes one multimodal file from the processed ophys NWB, the stimulus table, the running, Facemap and Lightning Pose outputs, and the Xenium tables.
- `make_ophys_nwb.py`: makes the ophys part from the per-plane processed files, for sessions where the AIND pipeline NWB has placeholder metadata or does not exist. The data are the same as in the AIND NWB.

## Citation and contact

If you use the data, cite the dandiset: DANDI:002015 (see the "Cite as" box on the [dandiset page](https://dandiarchive.org/dandiset/002015)).

Contact: Omid Zobeiri and Anton Arkhipov, Allen Institute. For questions about the code, open an issue in this repository.

Funding: NIH 1 U01 MH 130907-02.
