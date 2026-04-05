#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>

#include "simulation_driver.hpp"

#include "common_types.hpp"
#include "structure_data.hpp"
#include "forcefields/activated_interaction_forcefield.hpp"
#include "potentials/cosine_angle_potential.hpp"
#include "reactions/association_simulator.hpp"
#include "misc/point_manipulations.hpp"
#include "misc/make_vector.hpp"


simulation_driver::simulation_driver(simulation_config const& config)
: _config(config)
, _store(config.sampling.output_filename)
, _random(config.sampling.random_seed)
, _structure({})
{
    setup();
    save_metadata();
}


void
simulation_driver::setup()
{
    _setup = {};

    for (auto const& chain_config : _config.chains) {
        _setup.chains.push_back({
            .start  = _setup.particle_count,
            .end    = _setup.particle_count + chain_config.length,
            .config = chain_config
        });
        _setup.particle_count += chain_config.length;
    }
    _setup.box = md::open_box { .particle_count = _setup.particle_count };
    _setup.container_sphere = { .center = {}, .radius = _config.environment.container_radius };

    _structure = structure_data(_setup.box);

    setup_particles();
    setup_association_simulator();
    setup_loop_extrusion_simulator();
    setup_loop_capture_simulator();
    setup_forcefield_pairwise();
    setup_forcefield_connectivity();
    setup_forcefield_associations();
    setup_forcefield_extruders();
    setup_forcefield_captures();
    setup_forcefield_container();
}


void
simulation_driver::setup_particles()
{
    for (std::size_t i = 0; i < _setup.particle_count; i++) {
        _system.add_particle({
            .mobility = _config.chain.monomer_mobility,
        });
    }
}


void
simulation_driver::setup_association_simulator()
{
    association_simulator associations({
        .site_count           = _system.particle_count(),
        .valency              = _config.association.valency,
        .association_distance = _config.association.association_distance,
        .association_rate     = _config.association.association_rate,
        .dissociation_rate    = _config.association.dissociation_rate,
    });

    for (auto const& chain : _setup.chains) {
        for (auto const& feature : chain.config.association_features) {
            for (std::size_t i = feature.site.start; i < feature.site.end; i++) {
                std::size_t const site = chain.start + i;

                if (auto val = feature.association) { associations.set_association_factor(site, *val); }
                if (auto val = feature.dissociation) { associations.set_dissociation_factor(site, *val); }
                if (auto val = feature.valency) { associations.set_valency(site, *val); }
            }
        }
    }

    _associations = std::make_shared<association_simulator>(associations);
    _structure.request_neighbor_list(_config.association.association_distance);
}


void
simulation_driver::setup_loop_extrusion_simulator()
{
    loop_extrusion_simulator extruders({
        .site_count       = _system.particle_count(),
        .loading_rate     = _config.loop_extrusion.loading_rate,
        .unloading_rate   = _config.loop_extrusion.unloading_rate,
        .extrusion_rate   = _config.loop_extrusion.extrusion_rate,
        .contraction_rate = _config.loop_extrusion.contraction_rate,
        .crossing_factor  = _config.loop_extrusion.crossing_factor,
        .max_distance     = _config.loop_extrusion.max_distance.value_or(INFINITY),
    });

    auto const foreach_site = [](auto const& feature, auto callback) {
        for (std::size_t const i : feature.site) {
            if (feature.direction & site_directions::left) {
                callback(i, loop_extrusion_simulator::directions::leftward);
            }
            if (feature.direction & site_directions::right) {
                callback(i, loop_extrusion_simulator::directions::rightward);
            }
        }
    };

    for (auto const& chain : _setup.chains) {
        for (auto const& feature : chain.config.loop_extrusion_features) {
            foreach_site(feature, [&](auto index, auto dir) {
                std::size_t const site = chain.start + index;
                if (auto val = feature.loading) { extruders.set_loading_factor(site, dir, *val); }
                if (auto val = feature.unloading) { extruders.set_unloading_factor(site, dir, *val); }
                if (auto val = feature.arrival) { extruders.set_arrival_factor(site, dir, *val); }
                if (auto val = feature.departure) { extruders.set_departure_factor(site, dir, *val); }
            });
        }

        // Prevent extruders in adjacent chains from entering this chain.
        //  >>>|         |<<<
        // ----|---------|----
        //     ^start     ^end
        extruders.set_arrival_factor(
            chain.start, loop_extrusion_simulator::directions::rightward, 0
        );
        extruders.set_arrival_factor(
            chain.end - 1, loop_extrusion_simulator::directions::leftward, 0
        );
    }

    _extruders = std::make_shared<loop_extrusion_simulator>(extruders);
}


void
simulation_driver::setup_loop_capture_simulator()
{
    loop_capture_simulator captures({
        .site_count       = _system.particle_count(),
        .loading_rate     = _config.loop_capture.loading_rate,
        .unloading_rate   = _config.loop_capture.unloading_rate,
        .diffusivity      = _config.loop_capture.diffusivity,
        .crossing_factor  = _config.loop_capture.crossing_factor,
        .capture_distance = _config.loop_capture.capture_distance,
        .capture_rate     = _config.loop_capture.capture_rate,
        .release_rate     = _config.loop_capture.release_rate,
        .traffic_rate     = _config.loop_capture.traffic_rate,
    });

    for (auto const& chain : _setup.chains) {
        for (auto const& feature : chain.config.loop_capture_features) {
            for (std::size_t i = feature.site.start; i < feature.site.end; i++) {
                std::size_t const site = chain.start + i;
                if (auto val = feature.loading) { captures.set_loading_factor(site, *val); }
                if (auto val = feature.unloading) { captures.set_unloading_factor(site, *val); }
            }
        }
        for (auto const& track : chain.config.loop_capture_tracks) {
            double const sign = (track.start < track.end ? +1 : -1);
            std::size_t const lower = std::min(track.start, track.end);
            std::size_t const upper = std::max(track.start, track.end);
            for (std::size_t i = lower; i <= upper; i++) {
                std::size_t const site = chain.start + i;
                captures.set_capture_factor(site, track.capture.value_or(1));
                captures.set_release_factor(site, track.release.value_or(1));
                captures.set_traffic_factor(site, track.traffic.value_or(1) * sign);
            }
        }
    }

    _captures = std::make_shared<loop_capture_simulator>(captures);
    _structure.request_neighbor_list(_config.loop_capture.capture_distance);
}


void
simulation_driver::setup_forcefield_pairwise()
{
    md::softcore_potential<2, 3> const repulsive_potential {
        .energy   = _config.chain.repulsive_energy,
        .diameter = _config.chain.repulsive_distance,
    };

    md::softcore_potential<8, 3> const attractive_potential {
        .energy   = -_config.chain.attractive_energy,
        .diameter = _config.chain.attractive_distance,
    };

    if (_config.chain.attractive_distance > 0 && _config.chain.attractive_energy > 0) {
        md::scalar const neighbor_distance = std::max(
            repulsive_potential.diameter, attractive_potential.diameter
        );

        _system.add_forcefield(
            md::make_neighbor_pairwise_forcefield(
                repulsive_potential + attractive_potential
            )
            .set_unit_cell(_setup.box)
            .set_neighbor_distance(neighbor_distance)
        );
    } else {
        _system.add_forcefield(
            md::make_neighbor_pairwise_forcefield(repulsive_potential)
            .set_unit_cell(_setup.box)
            .set_neighbor_distance(repulsive_potential.diameter)
        );
    }
}


void
simulation_driver::setup_forcefield_connectivity()
{
    auto bonds = _system.add_forcefield(
        md::make_bonded_pairwise_forcefield(
            [&, this](md::index, md::index) {
                return md::spring_potential {
                    .spring_constant      = _config.chain.spring_constant * _spring_factor,
                    .equilibrium_distance = _config.chain.spring_length,
                };
            }
        )
    );

    auto angles = _system.add_forcefield(
        md::make_bonded_triplewise_forcefield(
            md::cosine_angle_potential {
                .bending_energy  = _config.chain.angle_energy,
                .preferred_angle = _config.chain.angle_preference,
            }
        )
    );

    for (auto const& chain : _setup.chains) {
        bonds->add_bonded_range(chain.start, chain.end);
        if (_config.chain.angle_energy > 0) {
            angles->add_bonded_range(chain.start, chain.end);
        }

        for (auto const& loop : chain.config.static_loops) {
            bonds->add_bonded_pair(chain.start + loop.pair.site_1, chain.start + loop.pair.site_2);
        }
    }
}


void
simulation_driver::setup_forcefield_associations()
{
    md::softwell_potential const potential {
        .energy         = _config.association.association_energy,
        .decay_distance = _config.association.association_distance,
    };

    _system.add_forcefield(
        make_activated_interaction_forcefield(_associations, potential, _setup.box)
    );
}


void
simulation_driver::setup_forcefield_extruders()
{
    auto const potential = [&, this](md::index, md::index) {
        return md::spring_potential {
            .spring_constant      = _config.loop_extrusion.spring_constant * _spring_factor,
            .equilibrium_distance = _config.loop_extrusion.spring_length,
        };
    };

    _system.add_forcefield(
        make_activated_interaction_forcefield(_extruders, potential, _setup.box)
    );
}


void
simulation_driver::setup_forcefield_captures()
{
    auto const potential = [&, this](md::index, md::index) {
        return md::spring_potential {
            .spring_constant      = _config.loop_capture.spring_constant * _spring_factor,
            .equilibrium_distance = _config.loop_capture.spring_length,
        };
    };

    _system.add_forcefield(
        make_activated_interaction_forcefield(_captures, potential, _setup.box)
    );
}


void
simulation_driver::setup_forcefield_container()
{
    _system.add_forcefield(
        md::make_sphere_outward_forcefield(
            md::harmonic_potential {
                .spring_constant = _config.environment.container_spring_constant
            }
        )
        .set_sphere(_setup.container_sphere)
    );
}


void
simulation_driver::run()
{
    run_initialization_particles();
    run_initialization_associations();
    run_initialization_extruders();
    run_initialization_captures();

    for (phase_config const& phase : _config.sampling.phases) {
        run_sampling(phase);
    }
}


void
simulation_driver::run_initialization_particles()
{
    auto const positions = _system.view_positions();

    for (auto const& chain : _setup.chains) {
        auto const chain_length = chain.end - chain.start;
        auto const unit_step = _config.initialization.initial_step_monomers;
        auto const path_length = (chain_length + unit_step - 1) / unit_step;
        auto const path_step =
            _config.chain.spring_length
            * double(unit_step)
            * _config.initialization.initial_step_scale;

        std::uniform_real_distribution<md::scalar> coord(
            -_setup.container_sphere.radius, _setup.container_sphere.radius
        );
        md::point centroid;
        do {
            md::vector const shift = {coord(_random), coord(_random), coord(_random)};
            centroid = _setup.container_sphere.center + shift;
        } while (md::distance(centroid, _setup.container_sphere.center) > _setup.container_sphere.radius);

        std::vector<md::point> const coarse_path = generate_random_walk(path_length, path_step, _random);
        std::vector<md::point> chain_path = interpolate_points(coarse_path, chain_length);
        move_centroid(chain_path, centroid);

        std::ranges::copy(
            chain_path,
            positions.subview(chain.start, chain.end - chain.start).begin()
        );
    }
}


void
simulation_driver::run_initialization_associations()
{
    // No explicit initialization; equilibrated in the relaxation phase.
}


void
simulation_driver::run_initialization_extruders()
{
    if (!_config.initialization.preload_cohesin) {
        return;
    }

    loop_extrusion_simulator::snapshot_type initial_state;
    std::size_t next_id = 0;

    for (auto const& chain : _setup.chains) {

        std::vector<double> affinity_mods(chain.end - chain.start, 1);
        for (auto const& feature : chain.config.loop_extrusion_features) {
            for (std::size_t i = feature.site.start; i < feature.site.end; i++) {
                if (auto val = feature.loading) { affinity_mods[i] *= *val; }
                if (auto val = feature.unloading) { affinity_mods[i] /= *val; }
            }
        }

        for (std::size_t i = 0; i < chain.end - chain.start; i++) {
            std::size_t const site = chain.start + i;

            double const affinity = affinity_mods[i]
                * _config.loop_extrusion.loading_rate
                / _config.loop_extrusion.unloading_rate;

            // Limit initial loading up to one for each site.
            std::poisson_distribution<int> count_distr(affinity);
            if (std::isinf(affinity) || count_distr(_random) > 0) {
                initial_state.loops.push_back({ .id = next_id++, .site_1 = site, .site_2 = site });
            }
        }
    }

    _extruders->load_state(initial_state);
}


void
simulation_driver::run_initialization_captures()
{
    if (!_config.initialization.preload_cohesin) {
        return;
    }

    loop_capture_simulator::snapshot_type initial_state;
    std::size_t next_id = 0;

    for (auto const& chain : _setup.chains) {

        std::vector<double> affinity_mods(chain.end - chain.start, 1);
        for (auto const& feature : chain.config.loop_capture_features) {
            for (std::size_t i = feature.site.start; i < feature.site.end; i++) {
                if (auto val = feature.loading) { affinity_mods[i] *= *val; }
                if (auto val = feature.unloading) { affinity_mods[i] /= *val; }
            }
        }

        for (std::size_t i = 0; i < chain.end - chain.start; i++) {
            std::size_t const site = chain.start + i;

            double const affinity = affinity_mods[i]
                * _config.loop_capture.loading_rate
                / _config.loop_capture.unloading_rate;

            // Limit initial loading up to one for each site.
            std::poisson_distribution<int> count_distr(affinity);
            if (std::isinf(affinity) || count_distr(_random) > 0) {
                initial_state.cohesins.push_back({
                    .id            = next_id++,
                    .loaded_site   = site,
                    .captured_site = std::nullopt,
                });
            }
        }
    }
    initial_state.metadata.next_id = next_id;

    _captures->load_state(initial_state);
}


void
simulation_driver::run_sampling(phase_config const& phase)
{
    md::step const steps = phase.steps;
    md::step const sampling_interval = phase.sampling_interval.value_or(_config.sampling.sampling_interval);
    md::step const logging_interval = phase.logging_interval.value_or(_config.sampling.logging_interval);
    md::scalar const timestep = phase.timestep.value_or(_config.sampling.timestep);

    _spring_factor = phase.spring_factor.value_or(1);

    auto const callback = [&, this](md::step step) {
        if (logging_interval > 0 && step % logging_interval == 0) {
            show_progress(phase.name, step);
        }

        if (sampling_interval > 0 && step % sampling_interval == 0) {
            check_sanity();
            save_sample(phase.name, step);
        }

        _structure.update(_system.view_positions());
        _associations->step(timestep, _structure, _random);
        _extruders->step(timestep, _structure, _random);
        _captures->step(timestep, _structure, _random);
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.environment.temperature,
        .timestep    = timestep,
        .steps       = steps,
        .seed        = _random(),
        .callback    = callback,
    });
}


void
simulation_driver::show_progress(std::string const& phase_name, md::step step)
{
    std::time_t const wallclock_time = std::time(nullptr);
    double const energy = _system.compute_energy();
    std::size_t const associations = _associations->active_pairs().size();
    std::size_t const extruders = _extruders->active_pairs().size();
    std::size_t const captures = _captures->active_pairs().size();

    std::clog
        << std::put_time(std::localtime(&wallclock_time), "%F %T")
        << " [" << phase_name << "]"
        << " step: " << step
        << " | energy: " << std::setprecision(4) << energy
        << " | assocs: " << associations
        << " | x-loops: " << extruders
        << " | c-loops: " << captures
        << std::endl;
}


void
simulation_driver::save_sample(std::string const& phase_name, md::step step)
{
    _store.save_snapshot(phase_name, {
        .step         = step,
        .positions    = make_vector(_system.view_positions()),
        .associations = _associations->dump_state(),
        .extruders    = _extruders->dump_state(),
        .captures     = _captures->dump_state(),
    });
}


void
simulation_driver::save_metadata()
{
    std::vector<simulation_store::metadata_record::chain_record> chains;
    for (auto const& chain : _setup.chains) {
        chains.push_back({
            .start = chain.start,
            .end   = chain.end,
        });
    }

    std::vector<simulation_store::metadata_record::particle_record> particles;
    for (std::size_t site_index = 0; site_index < _associations->site_count(); site_index++) {
        auto const& site_data = _associations->get_site_data(site_index);
        particles.push_back({
            .valency = site_data.valency,
            .association_factor = site_data.association_factor,
            .dissociation_factor = site_data.dissociation_factor,
        });
    }

    _store.save_metadata({
        .config    = _config,
        .chains    = chains,
        .particles = particles,
    });
}


void
simulation_driver::check_sanity() const
{
    for (md::point const& position : _system.view_positions()) {
        if (!std::isfinite(position.vector().norm())) {
            throw std::runtime_error("Numerical instability");
        }
    }
}
