"""
Find and download the files of DANDI:002015.

The dandiset is embargoed until 2027-12-28. Until then you must have access to
the dandiset, and you must set your DANDI API key in the environment:

    export DANDI_API_KEY=<your key>     # from https://dandiarchive.org (click your initials)

File names:

    sub-<mouse>/sub-<mouse>_ses-<YYYYMMDDTHHMMSS>_behavior+image+ophys.nwb.zarr   (natural-movie session)
    sub-<mouse>/sub-<mouse>_ses-<YYYYMMDDTHHMMSS>_behavior+ophys.nwb.zarr         (drifting-grating sessions)
"""
import os
import re
from pathlib import Path

import pandas as pd

DANDISET_ID = "002015"
VERSION = "draft"
_NAME = re.compile(r"sub-(?P<subject>\d+)_ses-(?P<session>\d{8}T\d{6})_(?P<suffix>[a-z+]+)\.nwb\.zarr$")


def _parse(path):
    m = _NAME.search(str(path))
    if m is None:
        return None
    return dict(subject=m["subject"], session=m["session"],
                session_type="natural_movies" if "+image" in m["suffix"] else "drifting_gratings")


def _client():
    from dandi.dandiapi import DandiAPIClient
    return DandiAPIClient.for_dandi_instance("dandi", token=os.environ.get("DANDI_API_KEY"))


def list_sessions(subject=None):
    """Return a table of the NWB files in the dandiset (one row per session).

    Columns: subject, session, session_type, path (in the dandiset), size_gb.
    """
    rows = []
    with _client() as client:
        ds = client.get_dandiset(DANDISET_ID, VERSION)
        prefix = f"sub-{subject}/" if subject else ""
        for asset in ds.get_assets_with_path_prefix(prefix):
            info = _parse(asset.path)
            if info:
                rows.append(dict(**info, path=asset.path, size_gb=round(asset.size / 1e9, 2)))
    df = pd.DataFrame(rows)
    return df.sort_values(["subject", "session"]).reset_index(drop=True) if len(df) else df


def download_sessions(paths, out_dir="data"):
    """Download NWB files from the dandiset to out_dir.

    paths: one path or a list of paths in the dandiset (the 'path' column of
    list_sessions()), or 'sub-<mouse>' to download all 4 sessions of a mouse.
    A file goes to out_dir/<file>; a folder goes to out_dir/sub-<mouse>/<file>.
    Files that are already downloaded are skipped.
    Returns the list of local paths.
    """
    from dandi.download import download

    if isinstance(paths, (str, Path)):
        paths = [paths]
    out_dir = Path(out_dir)
    wanted = [str(q).strip("/") for q in paths]
    # a folder (sub-<mouse>) must end with '/'
    urls = [f"dandi://dandi/{DANDISET_ID}@{VERSION}/{q}" + ("" if q.endswith(".nwb.zarr") else "/")
            for q in wanted]
    download(urls, out_dir, existing="skip")
    names = [Path(q).name for q in wanted]
    return [p for p in local_sessions(out_dir).path
            if any(p.name == n or p.name.startswith(n + "_") for n in names)]


def local_sessions(root="data"):
    """Return a table of the NWB files under root (the same columns as list_sessions, with local paths)."""
    rows = [dict(**_parse(p), path=p) for p in sorted(Path(root).rglob("*.nwb.zarr")) if _parse(p)]
    return pd.DataFrame(rows, columns=["subject", "session", "session_type", "path"])
