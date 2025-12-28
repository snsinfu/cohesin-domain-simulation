#include <array>
#include <string>
#include <vector>

#include <jsoncons/json.hpp>

#include "simulation_config.hpp"


JSONCONS_N_MEMBER_TRAITS(
    phase_config,

    // Required fields
    2,
    name,
    steps,

    // Optional fields
    timestep,
    logging_interval,
    sampling_interval
)

JSONCONS_N_MEMBER_TRAITS(
    sampling_config,

    // Required fields
    0,

    // Optional fields
    timestep,
    logging_interval,
    sampling_interval,
    random_seed,
    output_filename,
    phases
)


JSONCONS_N_MEMBER_TRAITS(
    environment_config,

    // Required fields
    0,

    // Optional fields
    temperature,
    container_radius,
    container_spring_constant
)


JSONCONS_N_MEMBER_TRAITS(
    initialization_config,

    // Required fields
    0,

    // Optional fields
    extruder_preloading,
    initial_step_monomers,
    initial_step_scale
)


JSONCONS_N_MEMBER_TRAITS(
    chain_type_config,

    // Required fields
    0,

    // Optional fields
    repulsive_distance,
    repulsive_energy,
    attractive_distance,
    attractive_energy,
    spring_length,
    spring_constant,
    angle_preference,
    angle_energy,
    monomer_mobility
)


JSONCONS_N_MEMBER_TRAITS(
    association_type_config,

    // Required fields
    0,

    // Optional fields
    valency,
    association_distance,
    association_rate,
    dissociation_rate,
    association_energy
)


JSONCONS_N_MEMBER_TRAITS(
    extruder_type_config,

    // Required fields
    0,

    // Optional fields
    loading_rate,
    unloading_rate,
    extrusion_rate,
    contraction_rate,
    crossing_factor,
    spring_length,
    spring_constant,
    max_distance
)


JSONCONS_N_MEMBER_TRAITS(
    loop_capture_type_config,

    // Required fields
    0,

    // Optional fields
    loading_rate,
    unloading_rate,
    capture_distance,
    capture_rate,
    release_rate,
    crossing_factor,
    linear_diffusivity,
    spring_length,
    spring_constant
)


// For convenience, we allow site_range to be an interval object {start, end}
// or a single site index in the config. The latter is interpreted as a single
// element interval {index, index + 1}.
template<class Json>
struct jsoncons::json_type_traits<Json, site_range>
{
    using allocator_type = typename Json::allocator_type;

    static bool is(const Json& j) noexcept
    {
        return j.is_integer() || (j.is_object() && j.contains("start") && j.contains("end"));
    }

    static site_range as(const Json& j)
    {
        site_range value;
        if (j.is_number()) {
            value.start = j.template as<md::index>();
            value.end = value.start + 1;
        } else {
            value.start = j.at("start").template as<md::index>();
            value.end = j.at("end").template as<md::index>();
        }
        return value;
    }

    static Json to_json(site_range const& value, allocator_type alloc = {})
    {
        Json j(jsoncons::json_object_arg_t(), jsoncons::semantic_tag::none, alloc);
        j.try_emplace("start", value.start);
        j.try_emplace("end", value.end);
        return j;
    }
};


// site_pair is serialized as a length-two array.
template<class Json>
struct jsoncons::json_type_traits<Json, site_pair>
{
    using allocator_type = typename Json::allocator_type;

    static bool is(const Json& j) noexcept
    {
        return j.is_array() && j.size() == 2;
    }

    static site_pair as(const Json& j)
    {
        return {
            .site_1 = j[0].template as<md::index>(),
            .site_2 = j[1].template as<md::index>(),
        };
    }

    static Json to_json(site_pair const& value, allocator_type alloc = {})
    {
        Json j(jsoncons::json_array_arg_t(), jsoncons::semantic_tag::none, alloc);
        j.push_back(value.site_1);
        j.push_back(value.site_2);
        return j;
    }
};


JSONCONS_N_MEMBER_TRAITS(
    association_feature_config,

    // Required fields
    1,
    site,

    // Optional fields
    association,
    dissociation,
    valency
)


JSONCONS_ENUM_TRAITS(
    site_directions,

    // Choices
    both,
    left,
    right
)


JSONCONS_N_MEMBER_TRAITS(
    extruder_feature_config,

    // Required fields
    1,
    site,

    // Optional fields
    direction,
    loading,
    unloading,
    arrival,
    departure
)


JSONCONS_N_MEMBER_TRAITS(
    loop_capture_feature_config,

    // Required fields
    1,
    site,

    // Optional fields
    loading,
    unloading,
    capture,
    release
)


JSONCONS_N_MEMBER_TRAITS(
    static_loop_config,

    // Required fields
    1,
    pair
)


JSONCONS_N_MEMBER_TRAITS(
    chain_config,

    // Required fields
    1,
    length,

    // Optional fields
    association_features,
    extruder_features,
    loop_capture_features,
    static_loops,
    attributes
)


JSONCONS_N_MEMBER_TRAITS(
    simulation_config,

    // Required fields
    0,

    // Optional fields
    sampling,
    environment,
    initialization,
    chain,
    association,
    extruder,
    loop_capture,
    chains
)


simulation_config
parse_simulation_config(std::string const& text)
{
    simulation_config config = jsoncons::decode_json<simulation_config>(text);
    config.config_text = text;
    return config;
}


std::string
format_simulation_config(simulation_config const& config)
{
    std::string text;
    jsoncons::encode_json(config, text);
    return text;
}
