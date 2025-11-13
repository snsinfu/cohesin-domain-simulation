#pragma once

#include <cstdint>
#include <utility>

#include "hive_list.hpp"
#include "pair_array.hpp"


// Set-like data structure for unordered pairs of integers below a set value.
// It supports quick insertion, erasure, lookup and iteration.
class pair_bitmap
{
    // Since we use an O(n^2)-sized bitmap for indexing, we cannot handle a
    // large number of elements anyway. Use smaller index types internally to
    // save space.
    using index_type  = std::uint16_t;
    using cursor_type = std::uint32_t;

    static constexpr cursor_type sentinel = cursor_type(-1);

public:
    using size_type = std::size_t;
    using pair_type = std::pair<index_type, index_type>;
    using iterator  = hive_list<pair_type>::const_iterator;

    pair_bitmap(size_type n)
    : _index{n, sentinel}
    {
    }

    // Returns the number of pairs that are set in the bitmap.
    size_type count() const noexcept { return _count; }

    // Forward range interface for iterating over set pairs.
    auto begin() const noexcept { return _pairs.begin(); }
    auto end()   const noexcept { return _pairs.end(); }

    // Unsets all pairs.
    void clear()
    {
        _index.fill(sentinel);
        _pairs.clear();
    }

    // Tests whether the given index pair is set. The order of the indices is
    // irrelevant, but the indices must not have the same value.
    bool test(size_type i, size_type j) const
    {
        return _index.at(i, j) != sentinel;
    }

    // Sets the pair (i,j).
    void set(size_type i, size_type j)
    {
        if (!test(i, j)) {
            auto const it = _pairs.insert({ index_type(i), index_type(j) });
            _index.at(i, j) = cursor_type(it.cursor());
            _count++;
        }
    }

    // Unsets the pair (i,j).
    void unset(size_type i, size_type j)
    {
        if (auto const cur = _index.at(i, j); cur != sentinel) {
            _pairs.erase(_pairs.locate(cursor_type(cur)));
            _index.at(i, j) = sentinel;
            _count--;
        }
    }

private:
    size_type               _count = 0;
    hive_list<pair_type>    _pairs;
    pair_array<cursor_type> _index;
};
