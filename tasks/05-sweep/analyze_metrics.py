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
    frames: slice = slice(None, None),
    n_jobs: int = 1,
):
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    LOG.info("Analyzing %d files (%d workers)", len(input_filenames), n_jobs)

    with multiprocessing.Pool(n_jobs) as pool:
        results = [
            result for result in pool.imap_unordered(
                do_compute_metrics,
                (dict(filename=filename, frames=frames) for filename in input_filenames),
            )
            if result is not None
        ]

    LOG.info("Obtained %d valid metrics", len(results))
    LOG.info("Writing to %s", output_filename)

    make_table(results).write_csv(output_filename)


@dataclasses.dataclass
class Metrics:
    filename: str
    config: dict
    sweep_data: dict
    rg: np.ndarray
    msd_alpha: float
    msd_alpha_a: float
    msd_alpha_b: float


def make_table(results: list[Metrics]) -> pl.DataFrame:
    records: list[dict] = []

    for r in results:
        records.append(
            {
                **r.sweep_data,
                "rg": r.rg,
                "msd_alpha": r.msd_alpha,
                "msd_alpha_a": r.msd_alpha_a,
                "msd_alpha_b": r.msd_alpha_b,
            }
        )

    return pl.DataFrame(records)


def do_compute_metrics(kwargs: dict) -> Metrics | None:
    try:
        return compute_metrics(**kwargs)
    except Exception:
        LOG.exception("Error analyzing file %s", kwargs["filename"])
        return None


def compute_metrics(filename: str, frames: slice) -> Metrics:
    LOG.info("Analyzing file %s", filename)

    with h5py.File(filename, "r") as store:
        config, config_source = load_config(store)
        phase_store = store["phases/production"]
        snapshots = load_snapshots(phase_store, frames)

    meta = config_source["@meta"]
    sweep_data = meta["sweep_data"]
    positions_samples = np.array([s.positions for s in snapshots])

    chain, = config["chains"]
    chain_length = chain["length"]

    default_valency = config["association"]["valency"]
    selector_a = np.full(chain_length, False)
    for feature in chain.get("association_features", []):
        is_a = False
        match feature:
            case {"valency": valency}:
                is_a = (valency < default_valency)

            case {"association": _} | {"dissociation": _}:
                assoc = feature.get("association", 1)
                dissoc = feature.get("dissociation", 1)
                is_a = (assoc < dissoc)

        if is_a:
            site = feature["site"]
            start = site["start"]
            end = site["end"]
            selector_a[start:end] = True

    selector_b = ~selector_a

    # Rg
    rg = compute_rg(positions_samples).mean()

    # One-point MSD
    msd_paths = subtract_centroid(positions_samples)
    msd_alpha, _ = compute_msd_params(msd_paths, lag=MSD_LAG)
    msd_alpha_a, _ = compute_msd_params(msd_paths[:, selector_a], lag=MSD_LAG)
    msd_alpha_b, _ = compute_msd_params(msd_paths[:, selector_b], lag=MSD_LAG)

    return Metrics(
        filename=filename,
        config=config,
        sweep_data=sweep_data,
        rg=rg,
        msd_alpha=msd_alpha,
        msd_alpha_a=msd_alpha_a,
        msd_alpha_b=msd_alpha_b,
    )


def load_config(store: h5py.Group) -> tuple[dict, dict]:
    config = json.loads(store["metadata/config"][()])
    config_source = json.loads(store["metadata/config_source"][()])
    return config, config_source


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
    parser.add_argument("--frames", type=str)
    parser.add_argument("--output", dest="output_filename", type=str)
    parser.add_argument("-j", dest="n_jobs", type=int)
    parser.add_argument("input_filenames", metavar="traj", type=str, nargs="*")
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
