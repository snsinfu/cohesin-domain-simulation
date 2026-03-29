import argparse
import copy
import itertools
import json
import os.path

import numpy as np

from make_config import resolve_config


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
    sweep_specs = meta["sweep_specs"]

    sweep_set = list(sweep_params(sweep_specs))

    master_random = np.random.Generator(np.random.PCG64(config_seed))
    randoms = master_random.spawn(len(sweep_set))

    for config_id, swept_params in enumerate(sweep_set):
        config_params_used = {**config_params, **swept_params}
        config, info = resolve_config(
            config_template,
            **config_params_used,
            random=randoms[config_id],
        )
        info["config_id"] = config_id
        config["@meta"]["config_info"] = info
        config["@meta"]["config_params"] = config_params_used

        output_filename = os.path.join(output_dirname, f"config-{config_id}.json")
        with open(output_filename, "w") as file:
            json.dump(config, file)


def sweep_params(specs: dict[str, list[any]]) -> dict[str, any]:
    for values in itertools.product(*specs.values()):
        params = dict(zip(specs.keys(), values))
        yield params


if __name__ == "__main__":
    main()
