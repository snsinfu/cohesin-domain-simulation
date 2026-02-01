import argparse
import copy
import dataclasses
import json

import numpy as np


def main(
    *,
    chain_length: int,
    enhancer_density: float,
    enhancer_size: float,
    gene_count: int,
    gene_size: float,
    gene_size_min: int,
    master_seed: int | None = None,
):
    random = np.random.Generator(np.random.PCG64(master_seed))
    enhancer_random, gene_random = random.spawn(2)

    enhancers = generate_enhancers(
        chain_length,
        enhancer_density,
        avg_size=enhancer_size,
        random=enhancer_random,
    )

    genes = generate_gene(
        chain_length,
        gene_count,
        avg_size=gene_size,
        min_size=gene_size_min,
        random=gene_random,
    )

    association_features = []
    for enhancer in enhancers:
        association_features.append({
            "site": {"start": enhancer.start, "end": enhancer.end},
            "valency": 1,
        })

    loop_capture_features = [
        {
            "site": {"start": 0, "end": chain_length},
            "loading": 0,
        }
    ]
    for enhancer in enhancers:
        loop_capture_features.append({
            "site": {"start": enhancer.start, "end": enhancer.end},
            "loading": 1,
        })

    loop_capture_tracks = []
    for gene in genes:
        loop_capture_tracks.append({"start": gene.tss, "end": gene.tes})

    chain = {
        "length": chain_length,
        "association_features": association_features,
        "loop_capture_features": loop_capture_features,
        "loop_capture_tracks": loop_capture_tracks,
    }
    print(json.dumps(chain))


# ----------------------------------------------------------------------------

@dataclasses.dataclass
class Enhancer:
    start: int
    end: int


def generate_enhancers(
    chain_length: int,
    density: float,
    *,
    density_tol: float = 0,
    avg_size: float = 1,
    random: np.random.Generator,
) -> list[Enhancer]:
    while True:
        pattern = generate_binary_pattern(
            phi=density,
            island_size=avg_size,
            length=chain_length,
            random=random,
        )

        actual_density = sum(pattern) / chain_length
        if abs(actual_density - density) > density_tol:
            continue

        break # accept

    return [Enhancer(s, e) for s, e in derive_active_intervals(pattern)]


# ----------------------------------------------------------------------------

@dataclasses.dataclass
class Gene:
    tss: int
    tes: int


def generate_gene(
    chain_length: int,
    count: int,
    *,
    avg_size: float,
    min_size: int,
    random: np.random.Generator,
) -> list[Gene]:
    while True:
        pattern = generate_binary_pattern(
            phi=(avg_size * count / chain_length),
            island_size=avg_size,
            length=chain_length,
            random=random,
        )
        intervals = derive_active_intervals(pattern)

        if len(intervals) != count:
            continue

        sizes = [(e - s) for s, e in intervals]
        if min(sizes) < min_size:
            continue

        break # accept

    genes = [Gene(tss=s, tes=(e - 1)) for s, e in intervals]
    reverses = random.binomial(1, 0.5, size=len(genes))

    for gene, reverse in zip(genes, reverses):
        if reverse:
            gene.tss, gene.tes = gene.tes, gene.tss

    return genes


# ----------------------------------------------------------------------------

def generate_binary_pattern(
    phi: float,
    length: int,
    *,
    island_size: float = 2,
    random: np.random.Generator,
) -> list[int]:
    if phi < 0 or phi > 1:
        raise Exception("phi must be within [0, 1]")
    if length <= 0:
        raise Exception("length must be >= 1")
    if island_size < 1:
        raise Exception("island_size must be >= 1")

    if phi == 1:
        return [1] * length

    if phi > 0.5:
        inverse_pattern = generate_binary_pattern(
            1 - phi,
            length,
            island_size=island_size,
            random=random,
        )
        return [1 - x for x in inverse_pattern]

    # Markov chain that produces a binary pattern with given avg "up" fraction.
    down_prob = 1 / island_size
    up_prob = phi / (1 - phi) * down_prob

    x = int(random.uniform(0, 1) < phi)

    pattern = [x]
    for step in range(1, length):
        if x == 0:
            if random.uniform(0, 1) < up_prob:
                x = 1
        else:
            if random.uniform(0, 1) < down_prob:
                x = 0
        pattern.append(x)

    return pattern


def derive_active_intervals(pattern: list[int]) -> list[tuple[int, int]]:
    intervals = []
    offset = 0

    while offset < len(pattern):
        try:
            start = pattern.index(1, offset)
        except ValueError:
            break
        try:
            end = pattern.index(0, start + 1)
        except ValueError:
            end = len(pattern)
        intervals.append((start, end))
        offset = end

    return intervals


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", "-s", dest="master_seed", type=int, default=None)
    parser.add_argument("--enhancer-density", type=float, default=0.2)
    parser.add_argument("--enhancer-size", type=int, default=2)
    parser.add_argument("--gene-count", type=int)
    parser.add_argument("--gene-size", type=int, default=50)
    parser.add_argument("--gene-size-min", type=int, default=10)
    parser.add_argument("chain_length", type=int)
    return vars(parser.parse_args())


if __name__ == "__main__":
    main(**parse_args())
