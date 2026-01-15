#include <algorithm>
#include <array>
#include <numeric>
#include <ranges>
#include <utility>
#include <vector>

#include "loop_capture_simulator.hpp"
#include "monte_carlo.hpp"

#include <iostream>


loop_capture_simulator::loop_capture_simulator(config_type const& config)
: _config(config)
, _sites(config.site_count)
{
}


std::size_t
loop_capture_simulator::site_count() const
{
    return _sites.size();
}


sized_iterator_range<loop_capture_simulator::active_iterator>
loop_capture_simulator::active_pairs() const
{
    return sized_iterator_range(
        active_iterator(_cohesins.begin(), _cohesins.end()),
        active_iterator(_cohesins.end(), _cohesins.end()),
        _capture_count
    );
}


void
loop_capture_simulator::set_loading_factor(std::size_t site, double factor)
{
    _sites[site].loading_factor = factor;
}


void
loop_capture_simulator::set_unloading_factor(std::size_t site, double factor)
{
    _sites[site].unloading_factor = factor;
}


void
loop_capture_simulator::set_capture_factor(std::size_t site, double factor)
{
    _sites[site].capture_factor = factor;
}


void
loop_capture_simulator::set_traffic_factor(std::size_t site, double factor)
{
    _sites[site].traffic_factor = factor;
}


void
loop_capture_simulator::set_release_factor(std::size_t site, double factor)
{
    _sites[site].release_factor = factor;
}



void
loop_capture_simulator::set_arrival_factor(std::size_t site, double factor)
{
    _sites[site].arrival_factor = factor;
}


void
loop_capture_simulator::set_departure_factor(std::size_t site, double factor)
{
    _sites[site].departure_factor = factor;
}


void
loop_capture_simulator::step(double dt, structure_data const& structure, random_engine& random)
{
    debug_check_invariant();

    step_unloading(dt, random);
    step_loading(dt, random);
    step_release(dt, random);
    step_capture(dt, structure, random);
    step_diffusion(dt, random);
    step_traffic(dt, random);
}


void
loop_capture_simulator::step_loading(double dt, random_engine& random)
{
    for (std::size_t site = 0; site < _sites.size(); site++) {
        double const affinity =
            _sites[site].loading_factor *
            small_pow(_config.crossing_factor, _sites[site].occupancy);
        double const effective_rate = _config.loading_rate * affinity;

        if (poisson_process(effective_rate, dt, random)) {
            _cohesins.push_back({
                .id            = _next_id++,
                .loaded_site   = site,
                .captured_site = std::nullopt,
            });
            _sites[site].occupancy += 1;
        }
    }
}


void
loop_capture_simulator::step_unloading(double dt, random_engine& random)
{
    for (std::size_t cohesin_index = 0; cohesin_index < _cohesins.size(); ) {
        // Make a copy as _cohesins entry will be removed in-place.
        cohesin_data const cohesin = _cohesins[cohesin_index];

        double const effective_rate =
            _config.unloading_rate *
            _sites[cohesin.loaded_site].unloading_factor;

        if (poisson_process(effective_rate, dt, random)) {
            _sites[cohesin.loaded_site].occupancy--;
            if (cohesin.captured_site) {
                _capture_count--;
                _sites[*cohesin.captured_site].occupancy--;
            }
            std::swap(_cohesins[cohesin_index], _cohesins.back());
            _cohesins.pop_back();
        } else {
            cohesin_index++;
        }
    }
}


void
loop_capture_simulator::step_capture(
    double dt,
    structure_data const& structure,
    random_engine& random
)
{
    md::array_view<md::point const> const positions = structure.positions();

    struct candidate_data
    {
        md::index site_index;
        double    capture_rate;
    };
    std::vector<candidate_data> candidates;

    // When crossing_factor is not 1, effective capture rate changes
    // depending on the order of capture events. Thus, scan cohesins
    // randomly to reduce the order bias.
    std::vector<std::size_t> cohesin_order(_cohesins.size());
    std::iota(cohesin_order.begin(), cohesin_order.end(), std::size_t(0));
    std::ranges::shuffle(cohesin_order, random);

    for (std::size_t const cohesin_index : cohesin_order) {
        cohesin_data& cohesin = _cohesins[cohesin_index];

        if (cohesin.captured_site) {
            continue;
        }

        md::point const position = positions[cohesin.loaded_site];

        // Collect all sites that can be captured and select a random one
        // using Gillespie algorithm.
        candidates.clear();
        double event_rate = 0;

        structure.query_neighbors(
            position,
            [&](md::index site_index) {
                if (md::distance(position, positions[site_index]) > _config.capture_distance) {
                    return;
                }
                if (site_index == cohesin.loaded_site) {
                    return;
                }

                double const affinity =
                    _sites[site_index].capture_factor *
                    small_pow(_config.crossing_factor, _sites[site_index].occupancy);
                double const effective_rate = _config.capture_rate * affinity;

                candidates.push_back({
                    .site_index   = site_index,
                    .capture_rate = effective_rate,
                });
                event_rate += effective_rate;
            }
        );

        if (poisson_process(event_rate, dt, random)) {
            double const event_value =
                std::uniform_real_distribution<double>(0, event_rate)(random);
            double rate_sum = 0;

            for (candidate_data const& candidate : candidates) {
                rate_sum += candidate.capture_rate;
                if (event_value < rate_sum) {
                    _capture_count++;
                    _sites[candidate.site_index].occupancy++;
                    cohesin.captured_site = candidate.site_index;
                    break;
                }
            }
        }
    }
}


void
loop_capture_simulator::step_release(double dt, random_engine& random)
{
    for (cohesin_data& cohesin : _cohesins) {
        if (!cohesin.captured_site) {
            continue;
        }

        std::size_t const site = *cohesin.captured_site;
        double const event_rate = _config.release_rate * _sites[site].release_factor;

        if (poisson_process(event_rate, dt, random)) {
            cohesin.captured_site = std::nullopt;
            _sites[site].occupancy--;
        }
    }
}


void
loop_capture_simulator::step_diffusion(double dt, random_engine& random)
{
    auto const compute_arrival_factor = [&](site_data const& site) {
        return site.arrival_factor * small_pow(_config.crossing_factor, site.occupancy);
    };

    auto const displace = [&](cohesin_data& cohesin, std::size_t dest_site) {
        assert(dest_site < _config.site_count);
        _sites[cohesin.loaded_site].occupancy--;
        _sites[dest_site].occupancy++;
        cohesin.loaded_site = dest_site;
    };

    for (cohesin_data& cohesin : _cohesins) {
        // We assume that cohesin is immobile when it's capturing another site.
        if (cohesin.captured_site) {
            continue;
        }

        double const departure = _sites[cohesin.loaded_site].departure_factor;

        double const minus_arrival =
            cohesin.loaded_site > 0 ?
                compute_arrival_factor(_sites[cohesin.loaded_site - 1]) : 0;
        double const plus_arrival =
            cohesin.loaded_site + 1 < _sites.size() ?
                compute_arrival_factor(_sites[cohesin.loaded_site + 1]) : 0;

        // Use Gillespie algorithm to choose which direction to move to.
        double const minus_rate = _config.diffusivity * departure * minus_arrival;
        double const plus_rate = _config.diffusivity * departure * plus_arrival;
        double const event_rate = minus_rate + plus_rate;

        if (poisson_process(event_rate, dt, random)) {
            std::uniform_real_distribution<double> event_value(0, event_rate);
            if (event_value(random) < minus_rate) {
                displace(cohesin, cohesin.loaded_site - 1);
            } else {
                displace(cohesin, cohesin.loaded_site + 1);
            }
        }
    }
}


void
loop_capture_simulator::step_traffic(double dt, random_engine& random)
{
    auto const displace = [&](cohesin_data& cohesin, std::size_t dest_site) {
        assert(dest_site < _config.site_count);
        std::size_t const curr_site = *cohesin.captured_site;
        cohesin.captured_site = dest_site;
        _sites[curr_site].occupancy--;
        _sites[dest_site].occupancy++;
    };

    for (cohesin_data& cohesin : _cohesins) {
        if (!cohesin.captured_site) {
            continue;
        }

        std::size_t const curr_site = *cohesin.captured_site;
        double const speed = _config.traffic_rate * _sites[curr_site].traffic_factor;

        if (curr_site == 0 && speed < 0) {
            continue;
        }
        if (curr_site + 1 == _config.site_count && speed > 0) {
            continue;
        }

        std::size_t const dest_site = speed > 0 ? curr_site + 1 : curr_site - 1;
        double const arrival = small_pow(_config.crossing_factor, _sites[dest_site].occupancy);
        double const event_rate = std::fabs(speed) * arrival;

        if (poisson_process(event_rate, dt, random)) {
            displace(cohesin, dest_site);
        }
    }
}


loop_capture_simulator::snapshot_type
loop_capture_simulator::dump_state() const
{
    snapshot_type state = {};

    for (cohesin_data const& cohesin : _cohesins) {
        state.cohesins.push_back({
            .id            = cohesin.id,
            .loaded_site   = cohesin.loaded_site,
            .captured_site = cohesin.captured_site,
        });
    }

    state.metadata = {
        .next_id = _next_id,
    };

    return state;
}


void
loop_capture_simulator::load_state(snapshot_type const& state)
{
    auto copy = *this;

    copy._capture_count = 0;
    copy._cohesins.clear();

    for (auto const& cohesin : state.cohesins) {
        copy._cohesins.push_back({
            .id            = cohesin.id,
            .loaded_site   = cohesin.loaded_site,
            .captured_site = cohesin.captured_site,
        });
    }

    for (site_data& site : copy._sites) {
        site.occupancy = 0;
    }
    for (cohesin_data const& cohesin : copy._cohesins) {
        copy._sites[cohesin.loaded_site].occupancy++;
        if (cohesin.captured_site) {
            copy._capture_count++;
            copy._sites[*cohesin.captured_site].occupancy++;
        }
    }

    copy._next_id = state.metadata.next_id;

    *this = std::move(copy);
}


void
loop_capture_simulator::debug_check_invariant() const
{
#ifndef NDEBUG
    std::size_t actual_capture_count = 0;
    for (cohesin_data const& cohesin : _cohesins) {
        if (cohesin.captured_site) {
            actual_capture_count++;
        }
    }
    assert(_capture_count == actual_capture_count);

    std::vector<std::size_t> actual_occupancy(_sites.size());
    for (cohesin_data const& cohesin : _cohesins) {
        actual_occupancy.at(cohesin.loaded_site)++;
        if (cohesin.captured_site) {
            actual_occupancy.at(*cohesin.captured_site)++;
        }
    }
    for (std::size_t site = 0; site < _sites.size(); site++) {
        assert(_sites.at(site).occupancy == actual_occupancy.at(site));
    }
#endif
}
