import argparse
import copy
import dataclasses
import json
import os.path

import numpy as np


def main(
    *,
    input_filename: str = "config_template.json",
    output_dirname: str = "_configs",
):
    with open(input_filename) as file:
        config_template = json.load(file)

    meta = config_template["@meta"]
    config_params = meta["config_params"]
    config_seed = meta["config_seed"]
    config_count = meta["config_count"]

    random = np.random.Generator(np.random.PCG64(config_seed))

    for config_id in range(config_count):
        config, info = resolve_config(
            config_template,
            **config_params,
            random=random,
        )
        config["@meta"]["config_info"] = info

        output_filename = os.path.join(output_dirname, f"config-{config_id}.json")
        with open(output_filename, "w") as file:
            json.dump(config, file)


def resolve_config(
    config: dict,
    *,
    chain_length: int,
    acetylation_level: float,
    acetylation_size: float,
    acetylation_feature: dict,
    static_loop_count: int,
    static_loop_size_min: int,
    static_loop_size_max: int | str | None = None,
    static_loop_crossable: bool,
    gene_count: int,
    gene_size: float,
    gene_size_min: int,
    random: np.random.Generator,
) -> tuple[dict, dict]:
    config = copy.deepcopy(config)
    info = {}

    acetyl_random, gene_random, loop_random = random.spawn(3)

    match static_loop_size_max:
        case str("random"):
            info["static_loop_size_max"] = static_loop_size_max = int(
                loop_random.integers(static_loop_size_min, chain_length)
            )
        case str(spec):
            raise Exception(f"Unrecognized keyword '{spec}' for static_loop_size_max")

    islands = generate_acetyl_islands(
        chain_length,
        acetylation_level,
        avg_size=acetylation_size,
        random=acetyl_random,
    )

    genes = generate_genes(
        chain_length,
        count=gene_count,
        avg_size=gene_size,
        min_size=gene_size_min,
        random=gene_random,
    )

    loops = generate_island_loops(
        islands,
        static_loop_count,
        min_size=static_loop_size_min,
        max_size=static_loop_size_max,
        allow_crossing=static_loop_crossable,
        random=loop_random,
    )

    # Modulate association parameters on acetylated islands.
    association_features = []
    for island in islands:
        association_features.append(
            {
                "site": {"start": island.start, "end": island.end},
                **acetylation_feature,
            }
        )

    # Restrict loading of loop-capture cohesin to acetylated islands.
    loop_capture_features = [
        {"site": {"start": 0, "end": chain_length}, "loading": 0},
    ]
    for island in islands:
        loop_capture_features.append(
            {
                "site": {"start": island.start, "end": island.end},
                "loading": 1,
            }
        )

    # Genes define loop-capture traffic.
    loop_capture_tracks = [
        {"start": gene.tss, "end": gene.tes} for gene in genes
    ]

    # Define static loops on the chain.
    static_loops = [
        {"pair": [loop.site_1, loop.site_2]} for loop in loops
    ]

    config["chains"] = [
        {
            "length": chain_length,
            "association_features": association_features,
            "loop_capture_features": loop_capture_features,
            "loop_capture_tracks": loop_capture_tracks,
            "static_loops": static_loops,
        }
    ]

    return config, info


# ----------------------------------------------------------------------------

@dataclasses.dataclass
class Island:
    start: int
    end: int


def generate_acetyl_islands(
    chain_length: int,
    density: float,
    *,
    density_tol: float = 0,
    avg_size: float = 1,
    random: np.random.Generator,
) -> list[Island]:
    while True:
        pattern = generate_binary_pattern(
            chain_length,
            phi=density,
            island_size=avg_size,
            random=random,
        )

        actual_density = sum(pattern) / chain_length
        if abs(actual_density - density) > density_tol:
            continue

        break # accept

    return [Island(s, e) for s, e in derive_active_intervals(pattern)]


# ----------------------------------------------------------------------------

@dataclasses.dataclass
class Gene:
    tss: int
    tes: int


def generate_genes(
    chain_length: int,
    count: int,
    *,
    avg_size: float,
    min_size: int,
    random: np.random.Generator,
) -> list[Gene]:
    while True:
        pattern = generate_binary_pattern(
            chain_length,
            phi=(avg_size * count / chain_length),
            island_size=avg_size,
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

@dataclasses.dataclass
class Loop:
    site_1: int
    site_2: int


def loop_size(loop: Loop) -> int:
    return loop.site_2 - loop.site_1 + 1


def generate_island_loops(
    islands: list[Island],
    count: int,
    *,
    min_size: int | None = None,
    max_size: int | None = None,
    allow_crossing: bool = False,
    random: np.random.Generator,
) -> list[Loop]:
    loops = []

    def _check_crossing(pair, loops):
        return any(
            ((pair.site_1 < loop.site_1 and pair.site_2 > loop.site_1) or
             (pair.site_1 < loop.site_2 and pair.site_2 > loop.site_2))
            for loop in loops
        )

    while len(loops) < count:
        island_1 = random.choice(islands)
        island_2 = random.choice(islands)
        i = int(random.integers(island_1.start, island_1.end))
        j = int(random.integers(island_2.start, island_2.end))
        pair = Loop(min(i, j), max(i, j))

        size = loop_size(pair)
        if min_size is not None and size < min_size:
            continue
        if max_size is not None and size > max_size:
            continue
        if pair in loops:
            continue
        if not allow_crossing and _check_crossing(pair, loops):
            continue

        loops.append(pair)

    loops.sort(key=lambda pair: pair.site_1 - pair.site_2)
    return loops


## ----------------------------------------------------------------------------

def generate_binary_pattern(
    length: int,
    phi: float,
    *,
    island_size: float = 2,
    random: np.random.Generator,
) -> list[int]:
    if length <= 0:
        raise Exception("length must be >= 1")
    if phi < 0 or phi > 1:
        raise Exception("phi must be within [0, 1]")
    if island_size < 1:
        raise Exception("island_size must be >= 1")

    if phi == 1:
        return [1] * length

    if phi > 0.5:
        inverse_phi = 1 - phi
        inverse_island_size = island_size * inverse_phi / phi
        inverse_pattern = generate_binary_pattern(
            length,
            phi=inverse_phi,
            island_size=inverse_island_size,
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


if __name__ == "__main__":
    main()
