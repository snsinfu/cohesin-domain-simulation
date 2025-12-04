#pragma once

#include <cstddef>
#include <iterator>
#include <limits>
#include <vector>

#include <md.hpp>

#include "../common_types.hpp"
#include "../structure_data.hpp"
#include "../misc/sized_iterator_range.hpp"


class loop_extrusion_simulator
{
public:
    struct loop_data
    {
        std::size_t id     = 0;
        std::size_t site_1 = 0;
        std::size_t site_2 = 0;
    };

    enum class directions
    {
        leftward,
        rightward,
    };

    struct config_type
    {
        std::size_t site_count;
        double      loading_rate      = 0;
        double      unloading_rate    = 0;
        double      extrusion_rate    = 0;
        double      contraction_rate  = 0;
        double      crossing_factor   = 0;
        double      max_distance      = std::numeric_limits<double>::infinity();
    };

    struct snapshot_type
    {
        struct loop_record
        {
            std::size_t id;
            std::size_t site_1;
            std::size_t site_2;
        };
        struct metadata_record
        {
            std::size_t next_id;
        };
        std::vector<loop_record> loops;
        metadata_record          metadata;
    };

    class active_iterator;

    explicit loop_extrusion_simulator(config_type const& config);

    std::size_t                           site_count() const;
    sized_iterator_range<active_iterator> active_pairs() const;

    void set_loading_factor(std::size_t site, directions dir, double factor);
    void set_unloading_factor(std::size_t site, directions dir, double factor);
    void set_arrival_factor(std::size_t site, directions dir, double factor);
    void set_departure_factor(std::size_t site, directions dir, double factor);

    void step(double dt, structure_data const& structure, random_engine& random);

    snapshot_type dump_state() const;
    void          load_state(snapshot_type const& state);

private:
    struct site_params
    {
        double loading_factor   = 1;
        double unloading_factor = 1;
        double arrival_factor   = 1;
        double departure_factor = 1;
    };

    site_params& get_site_params(std::size_t site, directions dir);

    void step_loading(double dt, random_engine& random);
    void step_unloading(double dt, random_engine& random);
    void step_sliding(double dt, structure_data const& structure, random_engine& random);

private:
    struct site_data
    {
        site_params leftward;
        site_params rightward;
        unsigned    occupancy = 0;
    };

    config_type            _config;
    std::vector<site_data> _sites;
    std::vector<loop_data> _loops;
    std::size_t            _next_id = 0;
};


class loop_extrusion_simulator::active_iterator
{
public:
    struct value_type
    {
        std::size_t i;
        std::size_t j;
    };
    using pointer           = void;
    using reference         = value_type;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;
    using iterator_concept  = std::forward_iterator_tag;
    using base_iterator     = std::vector<loop_data>::const_iterator;

    active_iterator() = default;

    explicit active_iterator(base_iterator iter)
    : _iter(iter)
    {
    }

    bool operator==(active_iterator const& other) const
    {
        return _iter == other._iter;
    }

    bool operator!=(active_iterator const& other) const
    {
        return !operator==(other);
    }

    reference operator*() const
    {
        return { _iter->site_1, _iter->site_2 };
    }

    active_iterator& operator++()
    {
        ++_iter;
        return *this;
    }

    active_iterator operator++(int)
    {
        auto copy = *this;
        operator++();
        return copy;
    }

private:
    base_iterator _iter;
};
