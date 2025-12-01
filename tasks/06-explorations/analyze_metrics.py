import argparse
import dataclasses

import h5py
import numpy as np
import scipy.stats


MSD_LAG = 10


def main(
    *,
    filename: str,
    phase: str = "production",
    frames: slice = slice(None, None),
):
    with h5py.File(filename, "r") as store:
        phase_store = store["phases"][phase]
        snapshots = load_snapshots(phase_store, frames)

    positions_samples = np.array([s.positions for s in snapshots])
    chain_length = positions_samples.shape[1]

    # Rg
    rg = compute_rg(positions_samples).mean()

    # One-point MSD
    msd_paths = subtract_centroid(positions_samples)
    msd_alpha, _ = compute_msd_params(msd_paths, lag=MSD_LAG)

    # Two-point MSD
    pair_separation = chain_length // 2
    pair_msd_deltas = np.swapaxes(
        [
            positions_samples[:, i, :] - positions_samples[:, i + pair_separation, :]
            for i in range(0, chain_length - pair_separation - 1)
        ],
        0, 1,
    )
    pair_msd_alpha, _ = compute_msd_params(pair_msd_deltas, lag=MSD_LAG)

    print(f"Rg: {rg:.4g}")
    print(f"1p-MSD exp: {msd_alpha:.4g}")
    print(f"2p-MSD exp: {pair_msd_alpha:.4g}")


@dataclasses.dataclass
class Snapshot:
    positions: np.ndarray       # (*,3) float
    associations: np.ndarray    # (*,2) int


def load_snapshots(store: h5py.Group, frames: slice) -> list[Snapshot]:
    snapshots = []
    for step in store[".steps"][frames]:
        sample = store[step]
        positions = sample["positions"][:]
        associations = sample["associations"]["pairs"][:]
        snapshots.append(Snapshot(positions=positions, associations=associations))
    return snapshots


def compute_rg(positions: np.ndarray) -> np.ndarray:
    return np.sqrt(np.var(positions, axis=-2).sum(axis=-1))


def subtract_centroid(paths: np.ndarray) -> np.ndarray:
    return paths - paths.mean(axis=1, keepdims=True)


def compute_msd_params(paths: np.ndarray, lag: int) -> tuple[float, float]:
    msd = collect_msd(paths, lag).mean(1)
    x = np.log(np.arange(1, len(msd)))
    y = np.log(msd[1:])
    linreg = scipy.stats.linregress(x, y)
    return linreg.slope, np.exp(linreg.intercept)


def collect_msd(path: np.ndarray, lag: int) -> np.ndarray:
    squared_dists = np.zeros((lag,) + path.shape[1:-1])
    for t in range(0, len(path) - lag):
        deltas = path[t:t + lag] - path[t]
        squared_dists += (deltas ** 2).sum(-1)
    return squared_dists / (len(path) - lag)


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", type=str)
    parser.add_argument("--frames", type=str)
    parser.add_argument("filename", type=str)
    args = vars(parser.parse_args())
    reparse_slice(args, "frames")
    return remove_none(args)


def reparse_slice(args: dict, key: str) -> None:
    if spec := args.get(key):
        start_spec, end_spec = spec.split("-", 1)
        start = int(start_spec) if start_spec else None
        end = int(end_spec) if end_spec else None
        args[key] = slice(start, end)


def remove_none(args: dict) -> dict:
    return {key: value for key, value in args.items() if value is not None}


if __name__ == "__main__":
    main(**parse_args())
