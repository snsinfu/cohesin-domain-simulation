#include <stdexcept>

#include "structure_data.hpp"


structure_data::structure_data(box_type const& box)
: _box(box)
{
}


void
structure_data::request_neighbor_list(md::scalar distance)
{
    if (distance > _neighbor_distance.value_or(0)) {
        _neighbor_distance = distance;
        _neighbor_searcher = md::neighbor_searcher<box_type>(_box, distance);
    }
}


structure_data::neighbor_list_type const&
structure_data::neighbor_pairs(md::scalar distance) const
{
    if (distance > _neighbor_distance.value_or(0)) {
        throw std::runtime_error("neighbor list is not defined for given distance");
    }
    return _neighbor_list;
}


void
structure_data::update(md::array_view<md::point const> positions)
{
    _positions.assign(positions.begin(), positions.end());

    if (_neighbor_distance) {
        _neighbor_searcher->set_points(positions);
        _neighbor_list.update(positions, *_neighbor_distance, _box);
    }
}
