import argparse
import json
import logging

from dataclasses import dataclass

import h5py
import gsd
import gsd.hoomd
import numpy as np


LOG = logging.getLogger(__name__)
GSD_SCHEMA = dict(application="", schema="hoomd", schema_version=(1, 0))
DIMENSION = 3


def main(
    *,
    input_filename: str,
    output_filename: str,
    phase: str = "production",
):
    with h5py.File(input_filename, "r") as store:
        phase_store = store[phase]
        mod = CombinedMod(ChromatinMod(), CohesinMod())

        with gsd.fl.open(output_filename, "w", **GSD_SCHEMA) as output:
            with gsd.hoomd.HOOMDTrajectory(output) as traj:
                dump_trajectory(phase_store, traj, mod)


# ------------------------------------------------------------------------
# Data structures for passing around snapshot data

@dataclass
class ParticlesData:
    type_ids: list[int]
    type_names: list[str]
    positions: np.ndarray


@dataclass
class BondsData:
    type_ids: list[int]
    type_names: list[str]
    pairs: list[tuple[int]]


def combine_types(*sources: list[any]) -> tuple[list[int], list[str]]:
    type_ids = list(sources[0].type_ids)
    type_names = list(sources[0].type_names)

    for source in sources[1:]:
        id_offset = len(type_names)
        type_ids.extend(i + id_offset for i in source.type_ids)
        type_names.extend(source.type_names)

    return type_ids, type_names


def combine_particles_data(*sources: list[ParticlesData]) -> ParticlesData:
    type_ids, type_names = combine_types(*sources)
    return ParticlesData(
        type_ids=type_ids,
        type_names=type_names,
        positions=np.concat([source.positions for source in sources]),
    )


def combine_bonds_data(*sources: list[ParticlesData]) -> BondsData:
    type_ids, type_names = combine_types(*sources)
    return BondsData(
        type_ids=type_ids,
        type_names=type_names,
        pairs=sum((source.pairs for source in sources), []),
    )


# ------------------------------------------------------------------------
# Classes for extracting snapshots from our custom trajectory file

class TopologyMod:
    def derive_particles(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> ParticlesData:
        return ParticlesData(
            type_ids=[],
            type_names=[],
            positions=np.zeros((0, DIMENSION)),
        )

    def derive_bonds(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> BondsData:
        return BondsData(
            type_ids=[],
            type_names=[],
            pairs=[],
        )


class CombinedMod(TopologyMod):
    def __init__(self, *mods: list[TopologyMod]):
        self._mods = mods

    def derive_particles(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> ParticlesData:
        sources = (mod.derive_particles(metadata, snapshot) for mod in self._mods)
        return combine_particles_data(*sources)

    def derive_bonds(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> BondsData:
        sources = (mod.derive_bonds(metadata, snapshot) for mod in self._mods)
        return combine_bonds_data(*sources)


class ChromatinMod(TopologyMod):
    PARTICLE_TYPE = "chromatin"
    BOND_TYPE = "chromatin"

    def derive_particles(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> ParticlesData:
        positions = snapshot["positions"][:]
        return ParticlesData(
            type_ids=([0] * len(positions)),
            type_names=[self.PARTICLE_TYPE],
            positions=positions,
        )

    def derive_bonds(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> BondsData:
        chain_ranges = metadata["chains/ranges"][:]
        pairs = []
        for start, end in chain_ranges:
            pairs.extend((i, i + 1) for i in range(start, end - 1))
        return BondsData(
            type_ids=([0] * len(pairs)),
            type_names=[self.BOND_TYPE],
            pairs=pairs,
        )


class CohesinMod(TopologyMod):
    PARTICLE_TYPE = "cohesin"
    BOND_TYPE = "cohesin"

    def derive_particles(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> ParticlesData:
        positions = snapshot["positions"][:]
        cohesin_bonds = [(i, j) for i, j in snapshot["extruders/sites"]]
        cohesin_positions = np.array([
            (positions[i] + positions[j]) / 2 for i, j in cohesin_bonds
        ])
        return ParticlesData(
            type_ids=([0] * len(cohesin_positions)),
            type_names=[self.PARTICLE_TYPE],
            positions=cohesin_positions,
        )

    def derive_bonds(
        self,
        metadata: h5py.Group,
        snapshot: h5py.Group,
    ) -> BondsData:
        cohesin_bonds = [(i, j) for i, j in snapshot["extruders/sites"]]
        return BondsData(
            pairs=cohesin_bonds,
            type_ids=([0] * len(cohesin_bonds)),
            type_names=[self.BOND_TYPE],
        )


def dump_trajectory(
    store: h5py.Group,
    traj: gsd.hoomd.HOOMDTrajectory,
    mod: TopologyMod,
):
    metadata = store["metadata"]
    box_shape = derive_box_shape(metadata)

    # Dump all snapshots.
    for step in store[".steps"]:
        snapshot = store[step]
        particles = mod.derive_particles(metadata, snapshot)
        bonds = mod.derive_bonds(metadata, snapshot)

        frame_data = FrameData(
            step=int(step),
            box_shape=box_shape,
            particle_types=particles.type_ids,
            particle_type_names=particles.type_names,
            particle_positions=particles.positions,
            particle_attributes={},
            bond_types=bonds.type_ids,
            bond_type_names=bonds.type_names,
            bond_pairs=bonds.pairs,
        )
        traj.append(make_hoomd_frame(frame_data))


def derive_box_shape(metadata: h5py.Group) -> tuple[float, float, float]:
    config = load_config(metadata)
    radius = config["environment"]["container_radius"]
    diameter = 2 * radius
    return (diameter, diameter, diameter)


def load_config(metadata: h5py.Group) -> dict:
    return json.loads(metadata["config"][()])


@dataclass
class FrameData:
    step: int
    box_shape: tuple[float, float, float]
    particle_type_names: list[str]
    particle_types: np.ndarray
    particle_positions: np.ndarray
    particle_attributes: dict[str, np.ndarray]
    bond_type_names: list[str]
    bond_pairs: np.ndarray
    bond_types: np.ndarray


def make_hoomd_frame(data: FrameData) -> gsd.hoomd.Frame:
    frame = gsd.hoomd.Frame()

    frame.configuration = gsd.hoomd.ConfigurationData()
    frame.configuration.box = (*data.box_shape, 0, 0, 0)
    frame.configuration.step = data.step

    frame.particles = gsd.hoomd.ParticleData()
    frame.particles.types = data.particle_type_names
    frame.particles.position = data.particle_positions
    frame.particles.typeid = data.particle_types
    frame.particles.N = len(frame.particles.position)

    frame.bonds = gsd.hoomd.BondData(2)
    frame.bonds.types = data.bond_type_names
    frame.bonds.group = data.bond_pairs
    frame.bonds.typeid = data.bond_types
    frame.bonds.N = len(frame.bonds.group)

    # These entries are unused but need to be set.
    frame.angles = gsd.hoomd.BondData(3)
    frame.dihedrals = gsd.hoomd.BondData(4)
    frame.impropers = gsd.hoomd.BondData(4)
    frame.pairs = gsd.hoomd.BondData(2)
    frame.constraints = gsd.hoomd.ConstraintData()

    frame.state = {}
    frame.log = {}

    # Custom attributes
    for key, values in data.particle_attributes.items():
        frame.log[f"particles/{key}"] = values

    frame.validate()

    return frame


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", type=str)
    parser.add_argument("input_filename", type=str)
    parser.add_argument("output_filename", type=str)
    return _remove_none(vars(parser.parse_args()))


def _remove_none(d: dict) -> dict:
    return {key: value for key, value in d.items() if value is not None}


if __name__ == "__main__":
    main(**parse_args())
