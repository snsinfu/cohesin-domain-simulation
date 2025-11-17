#pragma once

#include <ranges>
#include <vector>


template<std::ranges::range R>
std::vector<std::ranges::range_value_t<R>>
make_vector(R const& r)
{
    using std::begin;
    using std::end;
    return std::vector<std::ranges::range_value_t<R>>(begin(r), end(r));
}
