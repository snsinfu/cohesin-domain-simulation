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

    meta = config_template["meta"]
    parameter_ranges = meta["parameter_ranges"]
    config_seed = meta["config_seed"]
    config_count = meta["config_count"]

    random.seed(config_seed)

    for config_id in range(config_count):
        params = generate_params(parameter_ranges)
        config = resolve_config(config_template, random_seed=config_id, **params)
        config["meta"]["config_data"] = {"config_id": config_id, **params}

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
    chain_length: int,
    association_valency: int,
    association_rate: float,
    association_constant: float,
    association_energy: float,
    extruder_loading_constant: float,
    random_seed: int,
) -> dict:
    config = copy.deepcopy(config)

    dissociation_rate = association_rate / association_constant
    unloading_rate = config["extruder"]["unloading_rate"]
    loading_rate = unloading_rate * extruder_loading_constant

    config["sampling"].update({
        "random_seed": random_seed,
    })

    config["association"].update({
        "valency": association_valency,
        "association_rate": association_rate,
        "dissociation_rate": dissociation_rate,
        "association_energy": association_energy,
    })

    config["extruder"].update({
        "loading_rate": loading_rate,
    })

    chain = {"length": chain_length}
    config["chains"] = [chain]

    return config


main()
