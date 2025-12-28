import argparse
import copy
import json
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
    config_seed = meta["config_seed"]
    config_count = meta["config_count"]
    config_params = meta["config_params"]

    random.seed(config_seed)

    for config_id in range(config_count):
        config = resolve_config(config_template, **config_params)
        config["@meta"]["config_data"] = {
            "config_id": config_id,
        }

        output_filename = f"{output_dirname}/config-{config_id}.json"
        with open(output_filename, "w") as output_file:
            json.dump(config, output_file)


def resolve_config(
    config: dict,
    *,
    domains: list[dict],
    loop_count: int,
    min_loop_size: int,
    max_loop_size: int,
) -> dict:
    config = copy.deepcopy(config)

    chain, = config["chains"]

    loops = []
    for domain in domains:
        start = domain["start"]
        end = domain["end"]

        for i in range(loop_count):
            loop_size = random.randint(min_loop_size, max_loop_size)
            site_1 = random.randint(start, end - loop_size - 1)
            site_2 = site_1 + loop_size + 1
            loops.append({"pair": [site_1, site_2]})

    chain["static_loops"] = loops

    return config


main()
