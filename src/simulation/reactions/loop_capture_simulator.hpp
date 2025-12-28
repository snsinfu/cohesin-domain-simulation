#pragma once

#include <cstddef>
#include <iterator>
#include <limits>
#include <vector>

#include <md.hpp>

#include "../common_types.hpp"
#include "../structure_data.hpp"
#include "../misc/sized_iterator_range.hpp"


class loop_capture_simulator
{
public:
    struct config_type
    {
        std::size_t site_count;
        double      loading_rate       = 0;
        double      unloading_rate     = 1;
        double      capture_distance   = 0;
        double      capture_rate       = 0;
        double      release_rate       = 1;
        double      crossing_factor    = 0;
        double      linear_diffusivity = 0;
    };

    struct snapshot_type
    {
        struct cohesin_record
        {
            std::size_t                id;
            std::size_t                loaded_site;
            std::optional<std::size_t> captured_site;
        };

        struct metadata_record
        {
            std::size_t next_id;
        };

        std::vector<cohesin_record> cohesins;
        metadata_record             metadata;
    };

    class active_iterator;

    explicit loop_capture_simulator(config_type const& config);

    std::size_t                           site_count() const;
    sized_iterator_range<active_iterator> active_pairs() const;

    void set_loading_factor(std::size_t site, double factor);
    void set_unloading_factor(std::size_t site, double factor);
    void set_capture_factor(std::size_t site, double factor);
    void set_release_factor(std::size_t site, double factor);
    void set_arrival_factor(std::size_t site, double factor);
    void set_departure_factor(std::size_t site, double factor);

    void step(double dt, structure_data const& structure, random_engine& random);

    snapshot_type dump_state() const;
    void          load_state(snapshot_type const& state);

private:
    struct site_data
    {
        unsigned occupancy        = 0;
        double   loading_factor   = 1;
        double   unloading_factor = 1;
        double   capture_factor   = 1;
        double   release_factor   = 1;
        double   arrival_factor   = 1;
        double   departure_factor = 1;
    };

    struct cohesin_data
    {
        std::size_t                id;
        std::size_t                loaded_site;
        std::optional<std::size_t> captured_site;
    };

    void step_unloading(double dt, random_engine& random);
    void step_loading(double dt, random_engine& random);
    void step_release(double dt, random_engine& random);
    void step_capture(double dt, structure_data const& structure, random_engine& random);
    void step_sliding(double dt, random_engine& random);

    void debug_check_invariant() const;

private:
    config_type               _config;
    std::vector<site_data>    _sites;
    std::vector<cohesin_data> _cohesins;
    std::size_t               _next_id = 0;
    std::size_t               _capture_count = 0;
};


class loop_capture_simulator::active_iterator
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
    using base_iterator     = std::vector<cohesin_data>::const_iterator;

    active_iterator() = default;

    explicit active_iterator(base_iterator iter, base_iterator end)
    : _iter(iter), _end(end)
    {
        skip_non_loops();
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
        assert(_iter->captured_site);
        return { _iter->loaded_site, *_iter->captured_site };
    }

    active_iterator& operator++()
    {
        ++_iter;
        skip_non_loops();
        return *this;
    }

    active_iterator operator++(int)
    {
        auto copy = *this;
        operator++();
        return copy;
    }

private:
    void skip_non_loops()
    {
        while (_iter != _end && !_iter->captured_site) {
            ++_iter;
        }
    }

private:
    base_iterator _iter;
    base_iterator _end;
};
