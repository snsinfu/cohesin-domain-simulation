import copy
import dataclasses
import json
import random


def main(
    *,
    input_filename: str = "config_template.json",
    output_dirname: str = "_configs",
):
    with open(input_filename, "r") as file:
        config_template = json.load(file)

    meta = config_template["@meta"]
    sweep_seed = meta["sweep_seed"]
    sweep_count = meta["sweep_count"]
    sweep_params = meta["sweep_params"]
    extra_params = meta["extra_params"]

    for config_id in range(sweep_count):
        config, resolved_data = resolve_config(
            config_template, **sweep_params, **extra_params,
        )
        config["@meta"]["sweep_data"] = {
            "config_id": config_id,
            **resolved_data,
        }

        output_filename = f"{output_dirname}/config-{config_id}.json"
        with open(output_filename, "w") as output_file:
            json.dump(config, output_file)


def resolve_config(
    config: dict,
    *,
    domain_size: int,
    domain_count: int,
    loop_boundary: bool,
    loop_count: int,
    loop_size_min: int,
    acetylation_island: int,
    acetylation_level: float,
    acetylation_error: float,
    acetylation_feature: dict,
    linker_length: int,
    linker_feature: dict,
) -> tuple[dict, dict]:
    resolved_data = {}

    domains = [
        generate_domain(
            domain_size,
            loop_boundary,
            loop_count,
            loop_size_min,
            acetylation_island,
            acetylation_level,
            acetylation_error,
        )
        for _ in range(domain_count)
    ]
    domain_chain = link_domains(domains, linker_length)

    resolved_data["domain_intervals"] = [
        (domain.start, domain.end) for domain in domain_chain.domains
    ]

    acetylation_levels = []
    for domain in domains:
        phi = sum(island.end - island.start for island in domain.islands) / domain.size
        acetylation_levels.append(phi)
    resolved_data["acetylation_levels"] = acetylation_levels

    # Define a chain based on the generated domain definitions.
    chain = {"length": domain_chain.length}

    chain["association_features"] = [
        {"site": {"start": island.start, "end": island.end}, **acetylation_feature}
        for domain in domain_chain.domains
        for island in domain.islands
    ]

    chain["association_features"] += [
        {"site": {"start": linker.start, "end": linker.end}, **linker_feature}
        for linker in domain_chain.linkers
    ]

    chain["static_loops"] = [
        {"pair": [loop.site_1, loop.site_2]}
        for domain in domain_chain.domains
        for loop in domain.loops
    ]

    config = copy.deepcopy(config)
    config["chains"] = [chain]
    return config, resolved_data


@dataclasses.dataclass
class Interval:
    start: int  # inclusive
    end: int    # non-inclusive


@dataclasses.dataclass
class Pair:
    site_1: int
    site_2: int


def loop_size(pair: Pair) -> int:
    return pair.site_2 - pair.site_1 + 1


@dataclasses.dataclass
class Domain:
    start: int
    end: int
    islands: list[Interval]
    loops: list[Pair]

    @property
    def size(self) -> int:
        return self.end - self.start


@dataclasses.dataclass
class DomainChain:
    length: int
    domains: list[Domain]
    linkers: list[Interval]


def generate_domain(
    domain_size: int,
    loop_boundary: bool,
    loop_count: int,
    loop_size_min: int,
    acetylation_island: int,
    acetylation_level: float,
    acetylation_error: float,
) -> Domain:
    while True:
        pattern = generate_acetylation_pattern(
            acetylation_level, domain_size, island_size=acetylation_island,
        )
        actual_level = compute_acetylation_level(pattern)
        if abs(actual_level - acetylation_level) <= acetylation_error:
            break
    islands = derive_acetylation_intervals(pattern)
    loops = generate_island_loops(
        islands,
        count=loop_count,
        min_size=loop_size_min,
        domain_size=domain_size,
        require_boundaries=loop_boundary,
    )
    return Domain(0, domain_size, islands, loops)


def generate_island_loops(
    islands: list[Interval],
    count: int,
    min_size: int,
    domain_size: int,
    require_boundaries: bool,
) -> list[Pair]:
    loops = []

    if require_boundaries:
        loops.extend(generate_boundary_loops(islands, min_size, domain_size))

    while len(loops) < count:
        island_1 = random.choice(islands)
        island_2 = random.choice(islands)
        i = random.randint(island_1.start, island_1.end - 1)
        j = random.randint(island_2.start, island_2.end - 1)
        pair = Pair(min(i, j), max(i, j))
        if loop_size(pair) < min_size:
            continue
        if pair in loops:
            continue
        loops.append(pair)

    loops.sort(key=lambda pair: -loop_size(pair))
    return loops


def generate_boundary_loops(
    islands: list[Interval],
    min_size: int,
    domain_size: int,
) -> list[Pair]:
    loops = []

    # Left end
    while True:
        paired_island = random.choice(islands + [Interval(domain_size - 1, domain_size)])
        i = 0
        j = random.randint(paired_island.start, paired_island.end - 1)
        pair = Pair(i, j)
        if loop_size(pair) >= min_size:
            break
    loops.append(pair)

    if pair.site_2 == domain_size - 1:
        return loops

    # Right end
    while True:
        paired_island = random.choice(islands)
        i = random.randint(paired_island.start, paired_island.end - 1)
        j = domain_size - 1
        pair = Pair(i, j)
        if loop_size(pair) >= min_size:
            break
    loops.append(pair)

    return loops


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

    # Markov chain for producing a binary pattern with given avg composition.
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


def derive_acetylation_intervals(pattern: list[int]) -> list[Interval]:
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
        intervals.append(Interval(start, end))
        offset = end

    return intervals


def compute_acetylation_level(pattern: list[int]) -> float:
    return sum(pattern) / len(pattern)


def link_domains(domains: list[Domain], linker_length: int) -> DomainChain:
    chain = DomainChain(length=0, domains=[], linkers=[])
    offset = 0

    for domain in domains:
        if offset > 0:
            chain.linkers.append(Interval(offset, offset + linker_length))
            offset += linker_length

        domain = copy.deepcopy(domain)
        domain.start += offset
        domain.end += offset

        for island in domain.islands:
            island.start += offset
            island.end += offset

        for loop in domain.loops:
            loop.site_1 += offset
            loop.site_2 += offset

        chain.domains.append(domain)
        offset += domain.size

    chain.length = offset

    return chain


if __name__ == "__main__":
    main()
