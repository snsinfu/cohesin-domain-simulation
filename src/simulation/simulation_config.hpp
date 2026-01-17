#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <md.hpp>

#include "misc/counting_iterator.hpp"


/** Sampling phase. */
struct phase_config
{
    std::string               name;
    md::step                  steps;
    std::optional<md::scalar> timestep;
    std::optional<md::step>   logging_interval;
    std::optional<md::step>   sampling_interval;
    std::optional<md::scalar> spring_factor;
};


/** Configuration object for the entire simulation. */
struct sampling_config
{
    md::scalar                timestep          = 1;
    md::step                  logging_interval  = 0;
    md::step                  sampling_interval = 0;
    std::uint64_t             random_seed       = 0;
    std::string               output_filename;
    std::vector<phase_config> phases;
};


/** Simulation container. */
struct environment_config
{
    md::scalar temperature               = 1;
    md::scalar container_radius          = 0;
    md::scalar container_spring_constant = 0;
};


/** Initialization */
struct initialization_config
{
    bool       preload_cohesin       = false;
    md::index  initial_step_monomers = 1;
    md::scalar initial_step_scale    = 1;
};


/** Polymer chain parameters. */
struct chain_type_config
{
    md::scalar repulsive_distance  = 1;
    md::scalar repulsive_energy    = 1;
    md::scalar attractive_distance = 0;
    md::scalar attractive_energy   = 0;
    md::scalar spring_length       = 1;
    md::scalar spring_constant     = 1;
    md::scalar angle_preference    = 0;
    md::scalar angle_energy        = 0;
    md::scalar monomer_mobility    = 1;
};


/** Bead-bead association parameters. */
struct association_type_config
{
    md::index  valency              = 0;
    md::scalar association_distance = 0;
    md::scalar association_rate     = 0;
    md::scalar dissociation_rate    = 0;
    md::scalar association_energy   = 0;
};


/** Loop extruder parameters. */
struct loop_extrusion_type_config
{
    md::scalar                loading_rate     = 0;
    md::scalar                unloading_rate   = 1;
    md::scalar                extrusion_rate   = 0;
    md::scalar                contraction_rate = 0;
    md::scalar                crossing_factor  = 1;
    md::scalar                spring_length    = 1;
    md::scalar                spring_constant  = 0;
    std::optional<md::scalar> max_distance;
};


/**
 * Loop capture parameters. During cohesin is loaded on a chain, the loaded
 * site is reversiblly associative against other sites on the chain.
 */
struct loop_capture_type_config
{
    md::scalar loading_rate       = 0;
    md::scalar unloading_rate     = 1;
    md::scalar diffusivity        = 0;
    md::scalar capture_distance   = 1;
    md::scalar capture_rate       = 0;
    md::scalar release_rate       = 1;
    md::scalar traffic_rate       = 0;
    md::scalar crossing_factor    = 1;
    md::scalar spring_length      = 0;
    md::scalar spring_constant    = 0;
};


/** Half-open interval for specifying a region on a chain. */
struct site_range
{
    md::index start = 0;
    md::index end   = 0;
};

inline counting_iterator<md::index> begin(site_range const& r) { return counting_iterator(r.start); }
inline counting_iterator<md::index> end(site_range const& r) { return counting_iterator(r.end); }


/** Index pair designating paired sites. */
struct site_pair
{
    md::index site_1 = 0;
    md::index site_2 = 0;
};


/** Chain feature for overriding association parameters. */
struct association_feature_config
{
    site_range                   site;
    std::optional<md::scalar>    association;
    std::optional<md::scalar>    dissociation;
    std::optional<std::uint32_t> valency;
};


/** Bitmasks for specifying the directionality of a feature on a chain. */
enum class site_directions : unsigned
{
    left  = 0b01,
    right = 0b10,
    both  = left | right,
};

inline bool operator&(site_directions value, site_directions mask) { return (unsigned(value) & unsigned(mask)) != 0; }

/** Chain feature for overriding extruder parameters. */
struct loop_extrusion_feature_config
{
    site_range                site;
    site_directions           direction = site_directions::both;
    std::optional<md::scalar> loading;
    std::optional<md::scalar> unloading;
    std::optional<md::scalar> arrival;
    std::optional<md::scalar> departure;
};


/**
 * Chain feature for overriding loop capture parameters. Models existence of
 * barriers for cohesin, CTCF, etc. on chromatin.
 */
struct loop_capture_feature_config
{
    site_range                site;
    std::optional<md::scalar> loading;
    std::optional<md::scalar> unloading;
    std::optional<md::scalar> arrival;
    std::optional<md::scalar> departure;
};


/**
 * Chain feature for defining transportation tracks for loop capturing cohesin.
 * Models genes. Cohesin captures the starting site (TSS), one-directionally
 * transported over the track, and stops at the ending site (TES).
 */
struct loop_capture_track_config
{
    md::index                 start;
    md::index                 end;
    std::optional<md::scalar> capture;
    std::optional<md::scalar> release;
    std::optional<md::scalar> traffic;
};


/**
 * Chain feature for defining additional bonds (static loops) on a chain.
 */
struct static_loop_config
{
    site_pair pair;
};


/** . */
struct chain_config
{
    md::index                                      length;
    std::vector<association_feature_config>        association_features;
    std::vector<loop_extrusion_feature_config>     loop_extrusion_features;
    std::vector<loop_capture_feature_config>       loop_capture_features;
    std::vector<loop_capture_track_config>         loop_capture_tracks;
    std::vector<static_loop_config>                static_loops;
    std::map<std::string, std::vector<md::scalar>> attributes;
};


/** . */
struct simulation_config
{
    sampling_config            sampling;
    environment_config         environment;
    initialization_config      initialization;
    chain_type_config          chain;
    association_type_config    association;
    loop_extrusion_type_config loop_extrusion;
    loop_capture_type_config   loop_capture;
    std::vector<chain_config>  chains;
    std::string                config_text;
};


/** Parses JSON representation of `simulation_config` structure. */
simulation_config parse_simulation_config(std::string const& text);

/** Formats `simulation_config` structure as a JSON string. */
std::string format_simulation_config(simulation_config const& config);
