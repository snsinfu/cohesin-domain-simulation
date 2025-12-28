import argparse
import dataclasses
import json
import logging
import multiprocessing

import h5py
import numpy as np
import polars as pl
import scipy.stats


LOG = logging.getLogger()
LOG_LEVEL = logging.INFO
LOG_FORMAT = "[%(asctime)s] %(message)s"
LOG_DATE_FORMAT = "%F %T"
MSD_LAG = 10


def main(
    *,
    input_filenames: list[str],
    output_filename: str,
    phase_key: str = "production",
    job_count: int = 1,
):
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    LOG.info("Analyzing %d files (%d workers)", len(input_filenames), job_count)

    with multiprocessing.Pool(job_count) as pool:
        results = [
            result for result in pool.imap_unordered(
                do_compute_metrics,
                (dict(filename=filename, phase_key=phase_key) for filename in input_filenames),
            )
            if result is not None
        ]

    LOG.info("Obtained %d valid metrics", len(results))
    LOG.info("Writing to %s", output_filename)

    make_table(results).write_csv(output_filename)


@dataclasses.dataclass
class Metrics:
    config_id: int
    msd_alpha: float
    separation_scores: np.ndarray


def make_table(results: list[Metrics]) -> pl.DataFrame:
    records: list[dict] = []

    for result in results:
        record = {
            "config_id": result.config_id,
            "msd_alpha": result.msd_alpha,
            "separation_score": result.separation_scores.mean(),
        }
        records.append(record)

    return pl.DataFrame(records)


def do_compute_metrics(kwargs: dict) -> Metrics | None:
    try:
        return compute_metrics(**kwargs)
    except Exception:
        LOG.exception("Error analyzing file %s", kwargs["filename"])
        return None


def compute_metrics(filename: str, phase_key: str) -> Metrics:
    LOG.info("Analyzing file %s", filename)

    with h5py.File(filename, "r") as store:
        config, config_source = load_config(store)
        phase_store = store["phases"][phase_key]
        snapshots = load_snapshots(phase_store, frames=slice(0, None))

    meta = config_source["@meta"]
    config_params = meta["config_params"]
    config_data = meta["config_data"]

    # FIXME
    msd_lag = 7

    positions_samples = np.array([s.positions for s in snapshots])
    msd_alpha, _ = compute_msd_params( subtract_centroid(positions_samples), lag=msd_lag)

    # FIXME
    separation_scores = []
    for positions in positions_samples:
        distance_matrix = compute_distance_matrix(positions)
        D_11 = mask_lower(distance_matrix[:200, :200], k=0).mean()
        D_12 = distance_matrix[:200, 200:].mean()
        D_22 = mask_lower(distance_matrix[200:, 200:], k=0).mean()
        separation_score = D_12 / (D_11 + D_22)
        separation_scores.append(separation_score)
    separation_scores = np.array(separation_scores)

    return Metrics(
        config_id=config_data["config_id"],
        msd_alpha=msd_alpha,
        separation_scores=separation_scores,
    )


def mask_lower(matrix: np.ndarray, k: int) -> np.ma.MaskedArray:
    mask = np.tril(np.ones_like(matrix, dtype=bool), k=-k)
    return np.ma.masked_where(mask, matrix)


def compute_distance_matrix(positions: np.ndarray) -> np.ndarray:
    return np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)


def load_config(store: h5py.Group) -> tuple[dict, dict]:
    metadata = store["metadata"]
    config = json.loads(metadata["config"][()])
    config_source = json.loads(metadata["config_source"][()])
    return config, config_source


@dataclasses.dataclass
class Snapshot:
    positions: np.ndarray


def load_snapshots(store: h5py.Group, frames: slice) -> list[Snapshot]:
    snapshots = []
    for step in store[".steps"][frames]:
        sample = store[step]
        positions = sample["positions"][:]
        snapshots.append(Snapshot(positions=positions))
    return snapshots


def subtract_centroid(paths: np.ndarray) -> np.ndarray:
    return paths - paths.mean(axis=1, keepdims=True)


def compute_msd_params(paths: np.ndarray, lag: int) -> tuple[float, float]:
    msd = compute_msd(paths, lag).mean(1)
    x = np.log(np.arange(1, len(msd)))
    y = np.log(msd[1:])
    linreg = scipy.stats.linregress(x, y)
    return linreg.slope, np.exp(linreg.intercept)


def compute_msd(path: np.ndarray, lag: int) -> np.ndarray:
    squared_dists = np.zeros((lag,) + path.shape[1:-1])
    for t in range(0, len(path) - lag):
        deltas = path[t:t + lag] - path[t]
        squared_dists += (deltas ** 2).sum(-1)
    return squared_dists / (len(path) - lag)


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", dest="phase_key", type=str)
    parser.add_argument("--output", dest="output_filename", type=str)
    parser.add_argument("--jobs", "-j", dest="job_count", type=int)
    parser.add_argument("input_filenames", metavar="traj", type=str, nargs="*")
    return remove_none(vars(parser.parse_args()))


def remove_none(args: dict) -> dict:
    return {key: value for key, value in args.items() if value is not None}


if __name__ == "__main__":
    main(**parse_args())
