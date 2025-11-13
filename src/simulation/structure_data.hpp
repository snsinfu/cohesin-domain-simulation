#pragma once

#include <optional>
#include <vector>

#include <md.hpp>


/**
 * Data structure for passing MD information to kinetic simulators.
 */
class structure_data
{
public:
    using box_type           = md::open_box;
    using neighbor_list_type = md::neighbor_list<box_type>;

    struct use_image_tag {};
    static inline constexpr use_image_tag const use_image = {};

    structure_data(box_type const& box);

    void                            update(md::array_view<md::point const> positions);
    void                            request_neighbor_list(md::scalar distance);
    md::index                       particle_count() const;
    md::array_view<md::point const> positions() const;
    md::point const&                position(md::index i) const;
    md::scalar                      distance(md::index i, md::index j) const;
    md::scalar                      distance(md::index i, md::index j, use_image_tag) const;
    md::scalar                      distance(md::index i, md::point const& point) const;
    md::scalar                      distance(md::index i, md::point const& point, use_image_tag) const;
    neighbor_list_type const&       neighbor_pairs(md::scalar distance) const;

    template<typename Callback>
    void query_neighbors(md::point const& point, Callback&& callback) const;

private:
    box_type                                       _box;
    std::optional<md::scalar>                      _neighbor_distance;
    std::vector<md::point>                         _positions;
    neighbor_list_type                             _neighbor_list;
    std::optional<md::neighbor_searcher<box_type>> _neighbor_searcher;
};


inline md::index
structure_data::particle_count() const
{
    return _positions.size();
}


inline md::array_view<md::point const>
structure_data::positions() const
{
    return _positions;
}


inline md::point const&
structure_data::position(md::index i) const
{
    return _positions[i];
}


inline md::scalar
structure_data::distance(md::index i, md::index j) const
{
    return md::distance(_positions[i], _positions[j]);
}


inline md::scalar
structure_data::distance(md::index i, md::index j, use_image_tag) const
{
    return _box.shortest_displacement(_positions[i], _positions[j]).norm();
}


inline md::scalar
structure_data::distance(md::index i, md::point const& point) const
{
    return md::distance(_positions[i], point);
}


inline md::scalar
structure_data::distance(md::index i, md::point const& point, use_image_tag) const
{
    return _box.shortest_displacement(_positions[i], point).norm();
}


template<typename Callback>
void
structure_data::query_neighbors(md::point const& point, Callback&& callback) const
{
    class callback_iterator
    {
    public:
        explicit callback_iterator(Callback& callback) : _callback(callback) {}
        callback_iterator& operator*() { return *this; }
        callback_iterator& operator++() { return *this; }
        callback_iterator& operator++(int) { return *this; }
        void operator=(md::index index) const { _callback(index); }
    private:
        Callback& _callback;
    };
    _neighbor_searcher->query(point, callback_iterator(callback));
}
