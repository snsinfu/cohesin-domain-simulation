#pragma once

#include <memory>
#include <string>
#include <vector>

#include <md.hpp>

#include "simulation_config.hpp"
#include "simulation_store.hpp"
#include "reactions/association_simulator.hpp"
#include "reactions/loop_extrusion_simulator.hpp"


class simulation_driver
{
public:
    explicit simulation_driver(simulation_config const& config);
    void     run();

private:
    struct simulation_setup
    {
        struct chain_setup
        {
            std::size_t  start;
            std::size_t  end;
            chain_config config;
        };
        std::size_t              particle_count;
        md::open_box             box;
        md::sphere               container_sphere;
        std::vector<chain_setup> chains;
    };

    void setup();
    void setup_particles();
    void setup_association_simulator();
    void setup_loop_extrusion_simulator();
    void setup_forcefield_pairwise();
    void setup_forcefield_connectivity();
    void setup_forcefield_associations();
    void setup_forcefield_extruders();
    void setup_forcefield_container();

    void run_initialization_particles();
    void run_initialization_associations();
    void run_initialization_extruders();
    void run_sampling(phase_config const& phase);

    void show_progress(std::string const& phase_name, md::step step);
    void save_sample(std::string const& phase_name, md::step step);
    void save_metadata();

private:
    simulation_config                         _config;
    simulation_setup                          _setup;
    simulation_store                          _store;
    md::system                                _system;
    random_engine                             _random;
    structure_data                            _structure;
    std::shared_ptr<association_simulator>    _associations;
    std::shared_ptr<loop_extrusion_simulator> _extruders;
};
