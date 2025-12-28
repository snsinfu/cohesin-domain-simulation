import argparse
import dataclasses
import json
import multiprocessing

import h5py
import numpy as np
import scipy.stats


MSD_LAG = 11


def main(
    *,
    trajectory_files: list[str],
    output_file: str,
    phase_key: str,
    worker_count: int,
    msd_lag: int,
):
    with multiprocessing.Pool(worker_count) as pool:
        results = [
            result for result in pool.imap_unordered(
                do_compute_metrics,
                (
                    dict(
                        filename=filename,
                        phase_key=phase_key,
                        msd_lag=msd_lag,
                    )
                    for filename in trajectory_files
                ),
            )
            if result is not None
        ]

    with h5py.File(output_file, "w") as store:
        store["config"] = json.dumps(results[0].config)
        store["site_msds"] = np.array([r.site_msds for r in results])
        store["separation_scores"] = np.array([r.separation_history.mean() for r in results])


@dataclasses.dataclass
class Metrics:
    filename: str
    config: dict
    site_msds: np.ndarray
    separation_history: np.ndarray


def do_compute_metrics(kwargs: dict) -> Metrics | None:
    try:
        return compute_metrics(**kwargs)
    except Exception:
        return None


def compute_metrics(filename: str, phase_key: str, msd_lag: int) -> Metrics:
    with h5py.File(filename, "r") as store:
        config = json.loads(store["metadata/config"][()])
        phase_store = store["phases"][phase_key]
        positions_samples = load_positions_samples(phase_store)

    site_msds = collect_msds(positions_samples, lag=msd_lag)

    # FIXME
    separation_history = []
    for positions in positions_samples:
        distance_matrix = compute_distance_matrix(positions)
        D_11 = mask_lower(distance_matrix[:200, :200], k=0).mean()
        D_12 = distance_matrix[:200, 200:].mean()
        D_22 = mask_lower(distance_matrix[200:, 200:], k=0).mean()
        separation = D_12 / (D_11 + D_22)
        separation_history.append(separation)
    separation_history = np.array(separation_history)

    return Metrics(
        filename=filename,
        config=config,
        site_msds=site_msds,
        separation_history=separation_history,
    )


def load_positions_samples(store: h5py.Group) -> np.ndarray:
    positions_samples = []
    for step_key in store[".steps"]:
        sample = store[step_key]
        positions_samples.append(sample["positions"][:])
    return np.array(positions_samples)


def collect_msds(
    paths: np.ndarray,  # (time, particle, dim)
    lag: int,
) -> np.ndarray:
    squared_dists = np.zeros((lag, paths.shape[1]))
    for t in range(0, paths.shape[0] - lag):
        delta = paths[t:t + lag] - paths[t]
        squared_dists += (delta ** 2).sum(-1)
    return squared_dists / (paths.shape[0] - lag)


def compute_distance_matrix(positions: np.ndarray) -> np.ndarray:
    return np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)


def mask_lower(matrix: np.ndarray, k: int) -> np.ma.MaskedArray:
    mask = np.tril(np.ones_like(matrix, dtype=bool), k=-k)
    return np.ma.masked_where(mask, matrix)


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-count", "-j", type=int, default=1)
    parser.add_argument("--output", "-o", dest="output_file", type=str)
    parser.add_argument("--phase-key", type=str, default="production")
    parser.add_argument("--msd-lag", type=int, default=MSD_LAG)
    parser.add_argument("trajectory_files", type=str, nargs="*")
    return vars(parser.parse_args())


if __name__ == "__main__":
    main(**parse_args())
