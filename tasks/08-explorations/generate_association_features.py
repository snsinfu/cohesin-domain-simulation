import argparse
import json
import random
import sys


def main(
    *,
    length: int,
    target_level: float,
    island_size: float,
    tolerance: float,
    feature_json: str | None,
    index_offset: int = 0,
    random_seed: int | None,
):
    random.seed(random_seed)

    feature = json.loads(feature_json) if feature_json else {}

    while True:
        pattern = generate_acetylation_pattern(target_level, length, island_size=island_size)
        intervals = derive_acetylation_intervals(pattern)
        actual_level = compute_acetylation_level(pattern)
        if abs(actual_level - target_level) <= tolerance:
            break

    pattern_text = "".join(map(str, pattern))
    sys.stderr.write(f"pattern: {pattern_text}\n")
    sys.stderr.write(f"target level: {target_level:g}\n")
    sys.stderr.write(f"actual level: {actual_level:g}\n")

    features = []
    for start, end in intervals:
        start += index_offset
        end += index_offset
        features.append({
            "site": (
                start if start + 1 == end else {"start": start, "end": end}
            ),
            **feature,
        })
    print(json.dumps(features))


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


def derive_acetylation_intervals(pattern: list[int]) -> list[slice]:
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


def compute_acetylation_level(pattern: list[int]) -> float:
    return sum(pattern) / len(pattern)


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tolerance", type=float, default=1)
    parser.add_argument("--island-size", type=float, default=2)
    parser.add_argument("--index-offset", type=int, default=0)
    parser.add_argument("--random-seed", type=int, default=None)
    parser.add_argument("--feature-json", type=str, default=None)
    parser.add_argument("length", type=int)
    parser.add_argument("target_level", type=float)
    return vars(parser.parse_args())


main(**parse_args())
