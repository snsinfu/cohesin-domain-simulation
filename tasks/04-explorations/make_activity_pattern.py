#!/usr/bin/env python3

import argparse
import random


def main(
    *,
    length: int,
    target_activity: float,
    island_size: float,
    index_offset: int = 0,
    random_seed: int | None,
):
    random.seed(random_seed)

    pattern = generate_activity_pattern(target_activity, length, island_size=island_size)
    ranges = derive_activity_ranges(pattern)
    actual_activity = compute_activity_level(pattern)

    print("".join(map(str, pattern)))
    print("ranges: " + " ".join(f"{r.start}-{r.stop}"for r in ranges))
    print(f"activity: {actual_activity:.4g}")


def generate_activity_pattern(
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
        inverse_pattern = generate_activity_pattern(1 - phi, length)
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


def derive_activity_ranges(pattern: list[int]) -> list[slice]:
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


def compute_activity_level(pattern: list[int]) -> float:
    return sum(pattern) / len(pattern)


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--island-size", type=float, default=2)
    parser.add_argument("--index-offset", type=int, default=0)
    parser.add_argument("--random-seed", type=int, default=None)
    parser.add_argument("length", type=int)
    parser.add_argument("target_activity", type=float)
    return vars(parser.parse_args())


main(**parse_args())
