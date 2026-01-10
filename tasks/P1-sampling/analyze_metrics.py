import argparse
import dataclasses
import json
import logging
import multiprocessing

import h5py
import numpy as np
import scipy.stats


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
    LOG.info("Writing to %s", output_file)

    with h5py.File(output_file, "w") as store:
        store["config"] = json.dumps(results[0].config)
        store["files"] = [r.filename for r in results]
        store["site_msds"] = [r.site_msds for r in results]
        store["domain_rgs"] = [r.domain_rgs for r in results]
        store["domain_distances"] = [r.domain_distances for r in results]
        store["separation_scores"] = [r.separation_scores for r in results]


@dataclasses.dataclass
class Metrics:
    filename: str
    config: dict
    site_msds: np.ndarray           # (time, particle)
    domain_rgs: np.ndarray          # (time, domain)
    domain_distances: np.ndarray    # (time, domain, domain)
    separation_scores: np.ndarray   # (time)


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
        phase_store = store["phases"][phase_key]
        positions_samples = load_positions_samples(phase_store)

    msd_paths = positions_samples
    if subtract_centroid:
        msd_paths = msd_paths - msd_paths.mean(axis=1, keepdims=True)
    site_msds = collect_msds(msd_paths, max_lag=max_msd_lag)

    # Domains
    sweep_data = config_source["@meta"]["sweep_data"]
    interval_1, interval_2 = sweep_data["domain_intervals"]
    slice_1 = slice(*interval_1)
    slice_2 = slice(*interval_2)

    domain_rgs = []
    domain_distances = []
    separation_scores = []

    for positions in positions_samples:
        rg_1 = compute_rg(positions[slice_1])
        rg_2 = compute_rg(positions[slice_2])
        domain_rgs.append((rg_1, rg_2))

        centroid_1 = positions[slice_1].mean(axis=0)
        centroid_2 = positions[slice_2].mean(axis=0)
        centroids = np.array([centroid_1, centroid_2])
        domain_distances.append(compute_distance_matrix(centroids))

        distance_matrix = compute_distance_matrix(positions)
        D_11 = mask_lower(distance_matrix[slice_1, slice_1], k=0).mean()
        D_12 = distance_matrix[slice_1, slice_2].mean()
        D_22 = mask_lower(distance_matrix[slice_2, slice_2], k=0).mean()
        separation = D_12 / (D_11 + D_22)
        separation_scores.append(separation)

    domain_rgs = np.array(domain_rgs)
    domain_distances = np.array(domain_distances)
    separation_scores = np.array(separation_scores)

    return Metrics(
        filename=filename,
        config=config_source,
        site_msds=site_msds,
        domain_rgs=domain_rgs,
        domain_distances=domain_distances,
        separation_scores=separation_scores,
    )


def load_positions_samples(store: h5py.Group) -> np.ndarray:
    positions_samples = []
    for step_key in store[".steps"]:
        sample = store[step_key]
        positions_samples.append(sample["positions"][:])
    return np.array(positions_samples)


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


def compute_distance_matrix(positions: np.ndarray) -> np.ndarray:
    return np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)


def mask_lower(matrix: np.ndarray, k: int) -> np.ma.MaskedArray:
    mask = np.tril(np.ones_like(matrix, dtype=bool), k=-k)
    return np.ma.masked_where(mask, matrix)


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-count", "-j", type=int, default=1)
    parser.add_argument("--output", "-o", dest="output_file", type=str)
    parser.add_argument("--phase", dest="phase_key", type=str, default=DEFAULT_PHASE)
    parser.add_argument("--msd-lag", dest="max_msd_lag", type=int, default=DEFAULT_MSD_LAG)
    parser.add_argument("--subtract-centroid", action="store_true", default=False)
    parser.add_argument("trajectory_files", type=str, nargs="*")
    return vars(parser.parse_args())


if __name__ == "__main__":
    main(**parse_args())
