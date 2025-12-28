import argparse
import json

import polars as pl


def main(*, config_files: list[str], output_file: str):
    sweep_table = []

    for filename in config_files:
        with open(filename, "r") as file:
            config = json.load(file)

        sweep_data = config["@meta"]["sweep_data"]
        sweep_table.append(sweep_data)

    pl.DataFrame(sweep_table).write_csv(output_file)


def parse_args() -> dict[str, any]:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", dest="output_file", type=str)
    parser.add_argument("config_files", type=str, nargs="*")
    return vars(parser.parse_args())


main(**parse_args())
