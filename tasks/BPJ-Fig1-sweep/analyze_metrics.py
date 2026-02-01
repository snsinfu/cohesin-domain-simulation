import argparse
import dataclasses
import json
import logging
import multiprocessing

import h5py
import numpy as np
import polars as pl

from lib.stats import fit_power_law


LOG = logging.getLogger(__name__)
LOG_LEVEL = logging.INFO
LOG_FORMAT = "[%(asctime)s] %(message)s"
LOG_DATE_FORMAT = "%F %T"

DEFAULT_PHASE = "production"
DEFAULT_MSD_LAG = 50    # Frames


def main(
    *,
    trajectory_files: list[str],
    output_file: str,
    dump_file: str,
    phase_key: str,
    worker_count: int,
    max_msd_lag: int,
    subtract_centroid: bool,
):
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    LOG.info("Analyzing %d files (%d workers)", len(trajectory_files), worker_count)

    with multiprocessing.Pool(worker_count) as pool:
        results = [
            result for result in pool.imap_unordered(
                do_compute_metrics,
                (
                    dict(
                        filename=filename,
                        phase_key=phase_key,
                        max_msd_lag=max_msd_lag,
                        subtract_centroid=subtract_centroid,
                    )
                    for filename in trajectory_files
                ),
            )
            if result is not None
        ]

    LOG.info("Obtained %d valid metrics", len(results))

    # Data dump
    LOG.info("Writing dump data to %s", dump_file)

    with h5py.File(dump_file, "w") as store:
        store["config_ids"] = [r.config_id for r in results]
        store["site_msds"] = [r.site_msds for r in results]
        store["mean_msds"] = [r.site_msds.mean(axis=1) for r in results]

    # Stats output
    LOG.info("Writing stats to %s", output_file)

    output_table = pl.DataFrame([
        {
            "config_id": r.config_id,
            "rg_mean": r.rg_mean,
            "rg_std": r.rg_std,
            "msd_exponent": r.msd_exponent,
            "msd_prefactor": r.msd_prefactor,
            **r.sweep_params,
        }
        for r in results
    ])
    output_table.write_csv(output_file)


@dataclasses.dataclass
class Metrics:
    filename: str
    config_id: int
    sweep_params: dict
    rg_mean: float
    rg_std: float
    site_msds: np.ndarray       # (lag, site)
    msd_exponent: float
    msd_prefactor: float


def do_compute_metrics(kwargs: dict) -> Metrics | None:
    try:
        return compute_metrics(**kwargs)
    except Exception:
        LOG.exception("Error analyzing file %s", kwargs["filename"])
        return None


def compute_metrics(
    filename: str,
    phase_key: str,
    max_msd_lag: int,
    subtract_centroid: bool,
) -> Metrics:
    LOG.info("Analyzing file %s", filename)

    with h5py.File(filename, "r") as store:
        config_used = json.loads(store["metadata/config"][()])
        config_source = json.loads(store["metadata/config_source"][()])
        trajectory = load_trajectory(store["phases"][phase_key])

    meta = config_source["@meta"]
    info = meta["config_info"]
    config_id = info["config_id"]
    sweep_params = meta["sweep_params"]

    sampling_config = config_used["sampling"]
    phase_configs = sampling_config["phases"]
    for phase_config in phase_configs:
        if phase_config["name"] == phase_key:
            phase_config = {**sampling_config, **phase_config}
            frame_duration = phase_config["timestep"] * phase_config["sampling_interval"]
            break

    # Rg
    chain_rgs = np.array([
        compute_rg(positions) for positions in trajectory.positions_history
    ])

    # MSD
    msd_paths = trajectory.positions_history
    if subtract_centroid:
        msd_paths = msd_paths - msd_paths.mean(axis=1, keepdims=True)
    site_msds = collect_msds(msd_paths, max_lag=max_msd_lag)

    y = site_msds.mean(axis=1)
    x = np.arange(len(y)) * frame_duration
    msd_alpha, msd_beta = fit_power_law(x, y)

    return Metrics(
        filename=filename,
        config_id=config_id,
        sweep_params=sweep_params,
        rg_mean=chain_rgs.mean(),
        rg_std=chain_rgs.std(),
        site_msds=site_msds,
        msd_exponent=msd_alpha,
        msd_prefactor=msd_beta,
    )


# ----------------------------------------------------------------------------

@dataclasses.dataclass
class LoopHistory:
    birth_frame: int
    sites_history: list[tuple[int, int]]


@dataclasses.dataclass
class Trajectory:
    frame_count: int
    positions_history: np.ndarray   # (frames, sites, 3)
    loop_extrusion_history: dict[int, LoopHistory]
    loop_capture_history: dict[int, LoopHistory]


def load_trajectory(store: h5py.Group) -> Trajectory:
    positions_history = []
    loop_extrusion_history = {}
    loop_capture_history = {}

    def _collect_loop_history(history, loop_store):
        for loop_id, (i, j) in zip(loop_store["ids"], loop_store["sites"]):
            if loop_id not in history:
                history[loop_id] = LoopHistory(
                    birth_frame=frame_index, sites_history=[],
                )
            history[loop_id].sites_history.append((i, j))

    for frame_index, step_key in enumerate(store[".steps"]):
        sample_store = store[step_key]
        positions_history.append(sample_store["positions"][:])
        _collect_loop_history(loop_extrusion_history, sample_store["extruders"])
        _collect_loop_history(loop_capture_history, sample_store["captures"])

    return Trajectory(
        frame_count=len(positions_history),
        positions_history=np.array(positions_history),
        loop_extrusion_history=loop_extrusion_history,
        loop_capture_history=loop_capture_history,
    )


# ----------------------------------------------------------------------------

def collect_msds(
    paths: np.ndarray,  # (time, particle, dim)
    max_lag: int,
) -> np.ndarray:
    frame_count = max_lag + 1
    squared_dists = np.zeros((frame_count, paths.shape[1]))
    for t in range(0, paths.shape[0] - frame_count):
        delta = paths[t:t + frame_count] - paths[t]
        squared_dists += (delta ** 2).sum(-1)
    return squared_dists / (paths.shape[0] - frame_count)


def compute_rg(positions: np.ndarray) -> np.ndarray:
    centroids = positions.mean(axis=-2, keepdims=True)
    return np.sqrt((positions - centroids).var(axis=-2).sum(axis=-1))


# ----------------------------------------------------------------------------

def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-count", "-j", type=int, default=1)
    parser.add_argument("--output", "-o", dest="output_file", type=str)
    parser.add_argument("--dump", "-d", dest="dump_file", type=str)
    parser.add_argument("--phase", dest="phase_key", type=str, default=DEFAULT_PHASE)
    parser.add_argument("--msd-lag", dest="max_msd_lag", type=int, default=DEFAULT_MSD_LAG)
    parser.add_argument("--subtract-centroid", action="store_true", default=False)
    parser.add_argument("trajectory_files", type=str, nargs="*")
    return vars(parser.parse_args())


if __name__ == "__main__":
    main(**parse_args())
