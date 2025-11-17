#pragma once

#include <string>
#include <vector>

#include <h5.hpp>
#include <md.hpp>

#include "simulation_config.hpp"
#include "reactions/association_simulator.hpp"
#include "reactions/loop_extrusion_simulator.hpp"


class simulation_store
{
public:
    struct metadata_record
    {
        struct chain_record
        {
            md::index start;
            md::index end;
        };

        simulation_config         config;
        std::vector<chain_record> chains;
    };

    struct snapshot_record
    {
        md::step                                step;
        std::vector<md::point>                  positions;
        association_simulator::snapshot_type    associations;
        loop_extrusion_simulator::snapshot_type extruders;
    };

    explicit simulation_store(std::string const& filename);

    void save_metadata(std::string const& phase_name, metadata_record const& metadata);
    void save_snapshot(std::string const& phase_name, snapshot_record const& snapshot);

private:
    void do_save_positions(std::string const& path_base, snapshot_record const& snapshot);
    void do_save_associations(std::string const& path_base, snapshot_record const& snapshot);
    void do_save_extruders(std::string const& path_base, snapshot_record const& snapshot);

private:
    h5::file _store;
};
