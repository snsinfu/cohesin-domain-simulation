import copy
import json
import os
import math
import random


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


def generate_params(spec: dict[str, list]) -> dict[str, any]:
    params = {}
    for key_spec, value_spec in spec.items():
        key, kind = key_spec.split(":")
        params[key] = generate_value(kind, value_spec)
    return params


def generate_value(kind: str, value_spec: list) -> any:
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
    acetylation_level: float,
    acetylation_size: int,
    acetylation_valency: int,
    unloading_rate: float,
    loading_constant: float,
) -> tuple[dict, dict]:
    config = copy.deepcopy(config)
    resolved_params = {}

    # Parameterize extruder
    loading_rate = unloading_rate * loading_constant

    config["extruder"].update({
        "loading_rate": loading_rate,
        "unloading_rate": unloading_rate,
    })

    # Define chain
    chain, = config["chains"]
    chain_length = chain["length"]

    acetylation_pattern = generate_activity_pattern(
        acetylation_level, chain_length, island_size=acetylation_size,
    )
    acetylation_islands = derive_activity_ranges(acetylation_pattern)

    association_features = []
    for island in acetylation_islands:
        association_features.append({
            "site": {"start": island.start, "end": island.stop},
            "valency": acetylation_valency,
        })

    chain["association_features"] = (
        association_features + chain["association_features"]
    )

    config["chains"] = [chain]
    resolved_params["acetylation_level"] = compute_activity_level(acetylation_pattern)

    return config, resolved_params



main()
