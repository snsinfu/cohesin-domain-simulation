#pragma once

#include <memory>
#include <type_traits>

#include <md.hpp>


/**
 * Forcefield for integrating pairwise interactions activated by an
 * independent simulation. See simulators in the reactions directory.
 */
template<typename Simulator, typename PotFun, typename Box>
class activated_interaction_forcefield : public md::forcefield
{
public:
    using simulator_type = Simulator;
    using potfun_type = PotFun;
    using box_type = Box;

    activated_interaction_forcefield(
        std::shared_ptr<simulator_type const> const& simulator,
        potfun_type const&                           potfun,
        box_type const&                              box
    );

    md::scalar compute_energy(md::system const& system) override;
    void       compute_force(md::system const& system, md::array_view<md::vector> forces) override;

private:
    std::shared_ptr<simulator_type const> _simulator;
    potfun_type                           _potfun;
    box_type                              _box;
};


template<typename S, typename P, typename B>
activated_interaction_forcefield(
    std::shared_ptr<S> const& simulator,
    P const&                  potfun,
    B const&                  box
)
-> activated_interaction_forcefield<std::decay_t<S>, P, B>;


template<typename S, typename P, typename B>
auto make_activated_interaction_forcefield(
    std::shared_ptr<S> const& simulator,
    P const&                  pot,
    B const&                  box
)
{
    return activated_interaction_forcefield(
        simulator,
        md::detail::make_pair_potential_factory(pot),
        box
    );
}


template<typename S, typename P, typename B>
activated_interaction_forcefield<S, P, B>::activated_interaction_forcefield(
    std::shared_ptr<simulator_type const> const& simulator,
    potfun_type const&                           potfun,
    box_type const&                              box
)
: _simulator{simulator}, _potfun{potfun}, _box{box}
{
}


template<typename S, typename P, typename B>
md::scalar
activated_interaction_forcefield<S, P, B>::compute_energy(
    md::system const& system
)
{
    auto const positions = system.view_positions();
    md::scalar energy = 0;

    for (auto const& pair : _simulator->active_pairs()) {
        auto const potential = _potfun(system, pair.i, pair.j);
        auto const r_ij = _box.shortest_displacement(positions[pair.i], positions[pair.j]);
        energy += potential.evaluate_energy(r_ij);
    }

    return energy;
}


template<typename S, typename P, typename B>
void
activated_interaction_forcefield<S, P, B>::compute_force(
    md::system const&          system,
    md::array_view<md::vector> forces
)
{
    auto const positions = system.view_positions();

    for (auto const& pair : _simulator->active_pairs()) {
        auto const potential = _potfun(system, pair.i, pair.j);
        auto const r_ij = _box.shortest_displacement(positions[pair.i], positions[pair.j]);
        auto const f_ij = potential.evaluate_force(r_ij);
        forces[pair.i] += f_ij;
        forces[pair.j] -= f_ij;
    }
}
