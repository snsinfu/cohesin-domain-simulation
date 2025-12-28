import copy
import json
import os
import random
import math


def main(
    *,
    input_filename: str = "config_template.json",
    output_dirname: str = "_configs",
):
    with open(input_filename, "r") as file:
        config_template = json.load(file)

    meta = config_template["@meta"]
    sweep_ranges = meta["sweep_ranges"]
    sweep_seed = meta["sweep_seed"]
    sweep_count = meta["sweep_count"]
    extra_params = meta["extra_params"]

    random.seed(sweep_seed)

    for config_id in range(sweep_count):
        params = generate_params(sweep_ranges)
        config, resolved_params = resolve_config(
            config_template, **params, **extra_params,
        )
        config["@meta"]["sweep_data"] = {
            "config_id": config_id,
            **params,
            **resolved_params,
        }

        output_filename = f"{output_dirname}/config-{config_id}.json"
        with open(output_filename, "w") as output_file:
            json.dump(config, output_file)


def generate_params(spec: dict[str, any]) -> dict[str, any]:
    params = {}
    for key_spec, value_spec in spec.items():
        key, kind = key_spec.split(":")
        params[key] = generate_value(kind, value_spec)
    return params


def generate_value(kind: str, value_spec: any) -> any:
    match kind, value_spec:
        case "int", [lower, upper]:
            return random.randint(lower, upper)

        case "real", [lower, upper]:
            return random.uniform(lower, upper)

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


def resolve_config(
    config: dict,
    *,
    chain_length: int,
    acetylation_level: float,
    acetylation_size: float,
    acetylation_valency: int,
    acetylation_capture_factor: float,
    loading_constant: float,
    unloading_rate: float,
    capture_constant: int,
    release_rate: float,
    crossing_factor: float,
    linear_diffusivity: int,
) -> tuple[dict, dict]:
    config = copy.deepcopy(config)

    # Parameterize loop capture reactions
    loading_rate = loading_constant * unloading_rate
    capture_rate = capture_constant * release_rate

    config["loop_capture"].update({
        "loading_rate":       loading_rate,
        "unloading_rate":     unloading_rate,
        "capture_rate":       capture_rate,
        "release_rate":       release_rate,
        "crossing_factor":    crossing_factor,
        "linear_diffusivity": linear_diffusivity,
    })

    # Define chain
    acetylation_pattern = generate_acetylation_pattern(
        acetylation_level, chain_length, island_size=acetylation_size,
    )
    acetylation_islands = derive_acetylation_ranges(acetylation_pattern)

    association_features = []
    loop_capture_features = []
    for island in acetylation_islands:
        site = {"start": island.start, "end": island.stop}
        association_features.append({
            "site": site, "valency": acetylation_valency,
        })
        loop_capture_features.append({
            "site": site, "capture": acetylation_capture_factor,
        })

    config["chains"] = [
        {
            "length": chain_length,
            "association_features": association_features,
            "loop_capture_features": loop_capture_features,
        }
    ]

    #
    resolved_params = {
        "acetylation_level": compute_acetylation_level(acetylation_pattern),
    }

    return config, resolved_params


def generate_acetylation_pattern(
    phi: float,
    length: int,
    *,
    island_size: float = 2,
) -> list[int]:
    if phi < 0 or phi >= 1:
        raise Exception("phi must be within [0, 1)")
    if length <= 0:
        raise Exception("length must be >= 1")
    if island_size < 1:
        raise Exception("island_size must be >= 1")

    if phi > 0.5:
        inverse_pattern = generate_acetylation_pattern(1 - phi, length)
        return [1 - x for x in inverse_pattern]

    # Run a Markov chain that should produce a pattern obeying phi in average.
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


def derive_acetylation_ranges(pattern: list[int]) -> list[slice]:
    ranges = []
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
        ranges.append(slice(start, end))
        offset = end

    return ranges


def compute_acetylation_level(pattern: list[int]) -> float:
    return sum(pattern) / len(pattern)


main()
