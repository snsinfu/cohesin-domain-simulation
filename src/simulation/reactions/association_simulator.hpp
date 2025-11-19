#pragma once

#include <cstddef>
#include <iterator>
#include <vector>

#include "../common_types.hpp"
#include "../structure_data.hpp"
#include "../misc/pair_bitmap.hpp"
#include "../misc/sized_iterator_range.hpp"


class association_simulator
{
public:
    class active_iterator;

    struct config_type
    {
        std::size_t site_count;
        std::size_t valency              = std::size_t(-1); // no limit
        double      association_distance = 0;
        double      association_rate     = 0;
        double      dissociation_rate    = 0;
    };

    struct snapshot_type
    {
        struct association_record
        {
            std::size_t site_1;
            std::size_t site_2;
        };
        std::vector<association_record> associations;
    };

    explicit association_simulator(config_type const& config);

    sized_iterator_range<active_iterator> active_pairs() const;

    void set_association_factor(std::size_t site, double factor);
    void set_dissociation_factor(std::size_t site, double factor);
    void set_valency(std::size_t site, std::size_t valency);

    void step(double dt, structure_data const& structure, random_engine& random);

    snapshot_type dump_state() const;
    void          load_state(snapshot_type const& state);
    void          clear_state();

private:
    void step_dissociation(double dt, structure_data const& structure, random_engine& random);
    void step_association(double dt, structure_data const& structure, random_engine& random);

    struct site_data
    {
        std::size_t valency             = std::size_t(-1); // no limit
        std::size_t occupancy           = 0;
        double      association_factor  = 1;
        double      dissociation_factor = 1;
    };
    using site_pair = std::pair<std::size_t, std::size_t>;

private:
    std::vector<site_data> _sites;
    pair_bitmap            _associations;
    config_type            _config;
    std::vector<site_pair> _candidates_buffer;
};


class association_simulator::active_iterator
{
    using base_iterator = pair_bitmap::iterator;

public:
    struct value_type
    {
        std::size_t i = 0;
        std::size_t j = 0;
    };
    using reference         = value_type;
    using pointer           = void;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    active_iterator() = default;

    explicit active_iterator(base_iterator base)
    : _base{base}
    {
    }

    bool operator==(active_iterator const& other) const
    {
        return _base == other._base;
    }

    bool operator!=(active_iterator const& other) const
    {
        return !operator==(other);
    }

    reference operator*() const
    {
        return {.i = _base->first, .j = _base->second};
    }

    active_iterator& operator++()
    {
        ++_base;
        return *this;
    }

    active_iterator operator++(int)
    {
        auto copy = *this;
        operator++();
        return copy;
    }

private:
    base_iterator _base;
};
