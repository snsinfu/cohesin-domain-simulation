#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>


// Data structure for associating a value of type T to each pair of integers
// below a set value. It is a "symmetric" two-dimensional array without
// diagonal elements.
template<typename T>
class pair_array
{
    // Note: The flattened index for the pair (i, j), where 0<=i<j<n, is
    //
    //   a(i, j) = sum(n-1-k, 0<=k<=i-1) + j - i + 1
    //           = i (n-1) - i(i-1)/2 + j - i - 1
    //           = i (2n-i-3) / 2 + j - 1
    //
    // and the total number of pairs is
    //
    //   a(n-2, n-1) + 1 = (n-1) (n-2) / 2 + n - 1 .
    //
    // There are no pairs if n = 0 or 1.

public:
    using value_type       = T;
    using reference        = T&;
    using const_reference  = T const&;
    using size_type        = std::size_t;

    pair_array() = default;

    // Creates a pair_array for pairs of indices less than n.
    //
    // Note that the array is empty if n is 0 or 1.
    explicit pair_array(size_type n, value_type init = {})
    : _n{n}
    {
        if (n > 1) {
            _store.resize((n - 1) * (n - 2) / 2 + n - 1, init);
        }
    }

    // Assign value to all pairs.
    void fill(value_type const& value)
    {
        std::fill(_store.begin(), _store.end(), value);
    }

    // Returns a reference to the element associated to the pair (i,j). The
    // order of the pair is irrelevant; at(i, j) returns the same reference
    // as at(j, i). Note that i must not be equal to j.
    reference at(size_type i, size_type j)
    {
        return _store[locate(i, j)];
    }

    const_reference at(size_type i, size_type j) const
    {
        return _store[locate(i, j)];
    }

private:
    size_type locate(size_type i, size_type j) const
    {
        assert(i < _n);
        assert(j < _n);

        if (i > j) {
            std::swap(i, j);
        }
        assert(i < j);

        return i * (2 * _n - i - 3) / 2 + j - 1;
    }

private:
    size_type      _n = 0;
    std::vector<T> _store;
};
