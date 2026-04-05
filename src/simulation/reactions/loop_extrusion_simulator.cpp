#include <algorithm>
#include <array>
#include <cassert>
#include <utility>

#include "loop_extrusion_simulator.hpp"
#include "monte_carlo.hpp"


loop_extrusion_simulator::loop_extrusion_simulator(config_type const& config)
: _config(config)
, _sites(config.site_count)
{
}


std::size_t
loop_extrusion_simulator::site_count() const
{
    return _sites.size();
}


sized_iterator_range<loop_extrusion_simulator::active_iterator>
loop_extrusion_simulator::active_pairs() const
{
    return sized_iterator_range(
        active_iterator(_loops.begin()),
        active_iterator(_loops.end()),
        _loops.size()
    );
}


void
loop_extrusion_simulator::set_loading_factor(std::size_t site, directions dir, double factor)
{
    get_site_params(site, dir).loading_factor = factor;
}


void
loop_extrusion_simulator::set_unloading_factor(std::size_t site, directions dir, double factor)
{
    get_site_params(site, dir).unloading_factor = factor;
}


void
loop_extrusion_simulator::set_arrival_factor(std::size_t site, directions dir, double factor)
{
    get_site_params(site, dir).arrival_factor = factor;
}


void
loop_extrusion_simulator::set_departure_factor(std::size_t site, directions dir, double factor)
{
    get_site_params(site, dir).departure_factor = factor;
}


loop_extrusion_simulator::site_params&
loop_extrusion_simulator::get_site_params(std::size_t site, directions dir)
{
    assert(site < _sites.size());
    switch (dir) {
    case directions::leftward:
        return _sites[site].leftward;
    case directions::rightward:
        return _sites[site].rightward;
    }
}


void
loop_extrusion_simulator::step(double dt, structure_data const& structure, random_engine& random)
{
    step_loading(dt, random);
    step_unloading(dt, random);
    step_sliding(dt, structure, random);
}


void
loop_extrusion_simulator::step_loading(double dt, random_engine& random)
{
    for (std::size_t i = 0; i < _sites.size(); i++) {
        auto const affinity = harmonic_rate(
            get_site_params(i, directions::leftward).loading_factor,
            get_site_params(i, directions::rightward).loading_factor
        );
        auto const crossing = small_pow(_config.crossing_factor, _sites[i].occupancy);
        auto const effective_rate = _config.loading_rate * affinity * crossing;
        if (poisson_process(effective_rate, dt, random)) {
            _loops.push_back({
                .id     = _next_id++,
                .site_1 = i,
                .site_2 = i,
            });
            _sites[i].occupancy += 2;
        }
    }
}


void
loop_extrusion_simulator::step_unloading(double dt, random_engine& random)
{
    for (std::size_t i = 0; i < _loops.size(); ) {
        auto const loop = _loops[i];
        auto const affinity = harmonic_rate(
            get_site_params(loop.site_1, directions::leftward).unloading_factor,
            get_site_params(loop.site_2, directions::rightward).unloading_factor
        );
        auto const effective_rate = _config.unloading_rate * affinity;
        if (poisson_process(effective_rate, dt, random)) {
            _sites[loop.site_1].occupancy--;
            _sites[loop.site_2].occupancy--;
            std::swap(_loops[i], _loops.back());
            _loops.pop_back();
        } else {
            i++;
        }
    }
}


void
loop_extrusion_simulator::step_sliding(
    double dt,
    structure_data const& structure,
    random_engine& random
)
{
    auto const compute_affinity_factor = [&](site_params const& depart, site_params const& land) {
        return depart.departure_factor * land.arrival_factor;
    };

    auto const compute_crossing_factor = [&](site_data const& site) {
        return small_pow(_config.crossing_factor, site.occupancy);
    };

    auto const attempt_move = [&](
        std::size_t peer, std::size_t& stem, std::size_t dest, double rate
    ) {
        if (structure.distance(peer, dest) > _config.max_distance) {
            return;
        }
        auto const crossing = compute_crossing_factor(_sites[dest]);
        auto const effective_rate = rate * crossing;
        if (poisson_process(effective_rate, dt, random)) {
            _sites[stem].occupancy--;
            _sites[dest].occupancy++;
            stem = dest;
        }
    };

    auto const extrude_minus = [&](loop_data& loop) {
        if (loop.site_1 == 0) {
            return;
        }
        auto const affinity = compute_affinity_factor(
            get_site_params(loop.site_1, directions::leftward),
            get_site_params(loop.site_1 - 1, directions::leftward)
        );
        attempt_move(loop.site_2, loop.site_1, loop.site_1 - 1, _config.extrusion_rate * affinity);
    };

    auto const extrude_plus = [&](loop_data& loop) {
        if (loop.site_2 + 1 == _sites.size()) {
            return;
        }
        auto const affinity = compute_affinity_factor(
            get_site_params(loop.site_2, directions::rightward),
            get_site_params(loop.site_2 + 1, directions::rightward)
        );
        attempt_move(loop.site_1, loop.site_2, loop.site_2 + 1, _config.extrusion_rate * affinity);
    };

    auto const contract_minus = [&](loop_data& loop) {
        if (loop.site_1 == loop.site_2) {
            return;
        }
        auto const affinity = compute_affinity_factor(
            get_site_params(loop.site_1, directions::leftward),
            get_site_params(loop.site_1 + 1, directions::leftward)
        );
        attempt_move(loop.site_2, loop.site_1, loop.site_1 + 1, _config.contraction_rate * affinity);
    };

    auto const contract_plus = [&](loop_data& loop) {
        if (loop.site_1 == loop.site_2) {
            return;
        }
        auto const affinity = compute_affinity_factor(
            get_site_params(loop.site_2, directions::rightward),
            get_site_params(loop.site_2 - 1, directions::rightward)
        );
        attempt_move(loop.site_1, loop.site_2, loop.site_2 - 1, _config.contraction_rate * affinity);
    };

    for (auto& loop : _loops) {
        // Extrusion/contraction on minus/plus stems. We need to randomize
        // the order of the actions to avoid artificial bias.
        enum class action_type { extrude_minus, extrude_plus, contract_minus, contract_plus };

        std::array<action_type, 4> actions = {
            action_type::extrude_minus,
            action_type::extrude_plus,
            action_type::contract_minus,
            action_type::contract_plus,
        };
        std::shuffle(actions.begin(), actions.end(), random);

        for (auto const action : actions) {
            switch (action) {
            case action_type::extrude_minus:
                extrude_minus(loop);
                break;
            case action_type::extrude_plus:
                extrude_plus(loop);
                break;
            case action_type::contract_minus:
                contract_minus(loop);
                break;
            case action_type::contract_plus:
                contract_plus(loop);
                break;
            }
        }
    }
}


loop_extrusion_simulator::snapshot_type
loop_extrusion_simulator::dump_state() const
{
    snapshot_type state = {};

    for (auto const& loop : _loops) {
        state.loops.push_back({
            .id     = loop.id,
            .site_1 = loop.site_1,
            .site_2 = loop.site_2,
        });
    }

    state.metadata = {
        .next_id = _next_id,
    };

    return state;
}


void
loop_extrusion_simulator::load_state(snapshot_type const& state)
{
    auto copy = *this;

    copy._loops.clear();

    for (auto const& loop : state.loops) {
        copy._loops.push_back({
            .id     = loop.id,
            .site_1 = loop.site_1,
            .site_2 = loop.site_2,
        });
    }

    for (auto& site : copy._sites) {
        site.occupancy = 0;
    }
    for (auto const& loop : copy._loops) {
        copy._sites[loop.site_1].occupancy++;
        copy._sites[loop.site_2].occupancy++;
    }

    copy._next_id = state.metadata.next_id;

    *this = std::move(copy);
}
