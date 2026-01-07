#include "simulation_store.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "misc/h5_buffer_traits.hpp"


simulation_store::simulation_store(std::string const& filename)
: _store(filename, "w")
{
}

void
simulation_store::save_metadata(metadata_record const& metadata)
{
    std::string const path_root = "metadata";
    std::string const path_config = path_root + "/config";
    std::string const path_config_source = path_root + "/config_source";
    std::string const path_chain_ranges = path_root + "/chains/ranges";

    _store.dataset<h5::str>(path_config).write(format_simulation_config(metadata.config));
    _store.dataset<h5::str>(path_config_source).write(metadata.config.config_text);

    std::vector<std::array<std::size_t, 2>> chain_ranges;
    for (auto const& chain : metadata.chains) {
        chain_ranges.push_back({chain.start, chain.end});
    }
    _store.dataset<h5::i32, 2>(path_chain_ranges).write(chain_ranges);
}

void
simulation_store::save_snapshot(std::string const& phase_name, snapshot_record const& snapshot)
{
    std::string const step_key = std::to_string(snapshot.step);
    std::string const path_base = "phases/" + phase_name + "/" + step_key;
    std::string const path_steps = "phases/" + phase_name + "/.steps";

    do_save_positions(path_base, snapshot);
    do_save_associations(path_base, snapshot);
    do_save_extruders(path_base, snapshot);
    do_save_captures(path_base, snapshot);

    // Update the steps array.
    std::vector<std::string> steps;
    auto steps_dataset = _store.dataset<h5::str, 1>(path_steps);
    if (steps_dataset) {
        steps_dataset.read_fit(steps);
    }
    steps.push_back(step_key);
    steps_dataset.write(steps);
}


void
simulation_store::do_save_positions(std::string const& path_base, snapshot_record const& snapshot)
{
    std::string const path_positions = path_base + "/positions";
    _store.dataset<h5::f32, 2>(path_positions).write(snapshot.positions);
}


void
simulation_store::do_save_associations(std::string const& path_base, snapshot_record const& snapshot)
{
    std::string const path_associations_pairs = path_base + "/associations/pairs";

    std::vector<std::array<std::size_t, 2>> pairs;
    for (auto const& association : snapshot.associations.associations) {
        pairs.push_back({ association.site_1, association.site_2 });
    }

    _store.dataset<h5::i32, 2>(path_associations_pairs).write(pairs);
}


void
simulation_store::do_save_extruders(std::string const& path_base, snapshot_record const& snapshot)
{
    std::string const path_extruders = path_base + "/extruders";
    std::string const path_extruders_ids = path_extruders + "/ids";
    std::string const path_extruders_sites = path_extruders + "/sites";

    std::vector<std::size_t> ids;
    std::vector<std::array<std::size_t, 2>> sites;
    for (auto const& loop : snapshot.extruders.loops) {
        ids.push_back(loop.id);
        sites.push_back({ loop.site_1, loop.site_2 });
    }

    _store.dataset<h5::i32, 1>(path_extruders_ids).write(ids);
    _store.dataset<h5::i32, 2>(path_extruders_sites).write(sites);
}


void
simulation_store::do_save_captures(std::string const& path_base, snapshot_record const& snapshot)
{
    std::string const path_captures = path_base + "/captures";
    std::string const path_captures_ids = path_captures + "/ids";
    std::string const path_captures_sites = path_captures + "/sites";

    std::vector<std::size_t> ids;
    std::vector<std::array<std::size_t, 2>> sites;
    for (auto const& cohesin : snapshot.captures.cohesins) {
        ids.push_back(cohesin.id);

        // Convention: the captured site entry is set equal to the loaded site
        // when the cohesin is not capturing any site.
        sites.push_back({
            cohesin.loaded_site,
            cohesin.captured_site.value_or(cohesin.loaded_site),
        });
    }

    _store.dataset<h5::i32, 1>(path_captures_ids).write(ids);
    _store.dataset<h5::i32, 2>(path_captures_sites).write(sites);
}
