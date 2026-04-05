#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

#include "association_simulator.hpp"
#include "monte_carlo.hpp"


namespace
{
    using site_data = association_simulator::site_data;

    inline double
    mix_association_factor(double factor_i, double factor_j) {
        return std::sqrt(factor_i * factor_j);
    }

    inline double
    mix_dissociation_factor(double factor_i, double factor_j) {
        return std::sqrt(factor_i * factor_j);
    }
}


association_simulator::association_simulator(config_type const& config)
: _sites(config.site_count, { .valency = config.valency })
, _associations(config.site_count)
, _config(config)
{
}


sized_iterator_range<association_simulator::active_iterator>
association_simulator::active_pairs() const
{
    return sized_iterator_range(
        active_iterator(_associations.begin()),
        active_iterator(_associations.end()),
        _associations.count()
    );
}


void
association_simulator::set_association_factor(std::size_t site, double factor)
{
    assert(site < _sites.size());
    _sites[site].association_factor = factor;
}


void
association_simulator::set_dissociation_factor(std::size_t site, double factor)
{
    assert(site < _sites.size());
    _sites[site].dissociation_factor = factor;
}


void
association_simulator::set_valency(std::size_t site, std::size_t valency)
{
    assert(site < _sites.size());
    _sites[site].valency = valency;
}


site_data const&
association_simulator::get_site_data(std::size_t site) const
{
    assert(site < _sites.size());
    return _sites[site];
}


std::size_t
association_simulator::site_count() const
{
    return _sites.size();
}



void
association_simulator::step(double dt, structure_data const& structure, random_engine& random)
{
    step_dissociation(dt, structure, random);
    step_association(dt, structure, random);
}


void
association_simulator::step_dissociation(double dt, structure_data const&, random_engine& random)
{
    for (auto const& [i, j] : _associations) {
        auto const factor = mix_dissociation_factor(
            _sites[i].dissociation_factor,
            _sites[j].dissociation_factor
        );
        auto const rate = _config.dissociation_rate * factor;

        if (poisson_process(rate, dt, random)) {
            _associations.unset(i, j);
            _sites[i].occupancy--;
            _sites[j].occupancy--;
        }
    }
}


void
association_simulator::step_association(double dt, structure_data const& structure, random_engine& random)
{
    if (_config.association_rate == 0 ||
        _config.association_distance == 0) {
        return; // Avoid the hot loop below if possible.
    }

    // We use a two-pass algorithm to select reactions to take under the
    // valency constraint. Store candidates in this vector, shuffle it,
    // then apply allowed associations.
    std::vector<site_pair> candidates = std::move(_candidates_buffer);
    candidates.clear();

    for (auto const& [i, j] : structure.neighbor_pairs(_config.association_distance)) {
        auto const factor = mix_association_factor(
            _sites[i].association_factor,
            _sites[j].association_factor
        );
        auto const rate = _config.association_rate * factor;
        auto const distance = structure.distance(i, j);

        if (distance <= _config.association_distance &&
            !_associations.test(i, j) &&
            poisson_process(rate, dt, random)) {
            candidates.push_back({i, j});
        }
    }

    std::ranges::shuffle(candidates, random);

    for (auto const& [i, j] : candidates) {
        if (_sites[i].occupancy == _sites[i].valency ||
            _sites[j].occupancy == _sites[j].valency) {
            continue;
        }
        _associations.set(i, j);
        _sites[i].occupancy++;
        _sites[j].occupancy++;
    }

    _candidates_buffer = std::move(candidates);
}


association_simulator::snapshot_type
association_simulator::dump_state() const
{
    snapshot_type state;
    for (auto const& [i, j] : _associations) {
        state.associations.push_back({ .site_1 = i, .site_2 = j });
    }
    return state;
}


void
association_simulator::load_state(snapshot_type const& state)
{
    clear_state();

    for (auto const& association : state.associations) {
        _associations.set(association.site_1, association.site_2);
        _sites[association.site_1].occupancy++;
        _sites[association.site_2].occupancy++;
    }
}


void
association_simulator::clear_state()
{
    for (site_data& site : _sites) {
        site.occupancy = 0;
    }
    _associations.clear();
}
