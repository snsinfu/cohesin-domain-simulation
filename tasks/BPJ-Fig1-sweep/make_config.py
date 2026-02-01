import argparse
import copy
import dataclasses
import json
import math
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
    sweep_ranges = meta["sweep_ranges"]
    config_params = meta["config_params"]
    config_seed = meta["config_seed"]
    config_count = meta["config_count"]

    master_random = np.random.Generator(np.random.PCG64(config_seed))
    randoms = master_random.spawn(config_count)

    for config_id in range(config_count):
        params_random, resolve_random = randoms[config_id].spawn(2)
        sweep_params = generate_params(sweep_ranges, random=params_random)
        config, info = resolve_config(
            config_template,
            **config_params,
            **sweep_params,
            random=resolve_random,
        )
        info["config_id"] = config_id
        config["@meta"]["config_info"] = info
        config["@meta"]["sweep_params"] = sweep_params

        output_filename = os.path.join(output_dirname, f"config-{config_id}.json")
        with open(output_filename, "w") as file:
            json.dump(config, file)


# ----------------------------------------------------------------------------

def generate_params(
    specs: dict[str, any],
    random: np.random.Generator,
) -> dict[str, any]:
    params = {}
    for key_spec, value_spec in specs.items():
        key, kind = key_spec.split(":")
        params[key] = generate_value(kind, value_spec, random=random)
    return params


def generate_value(
    kind: str,
    value_spec: any,
    random: np.random.Generator,
) -> any:
    match kind, value_spec:
        case "int", [lower, upper]:
            return int(random.integers(lower, upper + 1))

        case "real", [lower, upper]:
            return float(random.uniform(lower, upper))

        case "logint", [lower, upper]:
            log_value = random.uniform(math.log(lower), math.log(upper + 1))
            return int(math.exp(log_value))

        case "logreal", [lower, upper]:
            log_value = random.uniform(math.log(lower), math.log(upper))
            return math.exp(log_value)

        case "choice", _:
            return random.choice(values)

        case _:
            raise Exception(f"unrecognized kind '{kind}' with spec '{value_spec}'")


# ----------------------------------------------------------------------------

def resolve_config(
    config: dict,
    *,
    # association parameters
    association_valency: int | None = None,
    association_energy: float | None = None,
    dissociation_rate: float | None = None,
    # chromatin spec
    chain_length: int,
    acetylation_level: float = 0,
    acetylation_level_tol: float = 0,
    acetylation_size: float = 1,
    acetylation_feature: dict = {},
    # RNG
    random: np.random.Generator,
) -> tuple[dict, dict]:
    config = copy.deepcopy(config)
    info = {}

    islands = generate_acetyl_islands(
        chain_length,
        acetylation_level,
        density_tol=acetylation_level_tol,
        avg_size=acetylation_size,
        random=random,
    )
    info["islands"] = [(island.start, island.end) for island in islands]
    info["acetylation_level"] = sum(island.end - island.start for island in islands) / chain_length

    # Modulate association parameters on acetylated islands.
    association_features = []
    for island in islands:
        association_features.append(
            {
                "site": {"start": island.start, "end": island.end},
                **acetylation_feature,
            }
        )

    # Define chain
    config["chains"] = [
        {
            "length": chain_length,
            "association_features": association_features,
        }
    ]

    # Association parameters
    def set_association_parameter(key: str, value: any):
        if value is not None:
            config["association"][key] = value

    set_association_parameter("valency", association_valency)
    set_association_parameter("association_energy", association_energy)
    set_association_parameter("dissociation_rate", dissociation_rate)

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
        inverse_pattern = generate_binary_pattern(
            length,
            phi=(1 - phi),
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


if __name__ == "__main__":
    main()
