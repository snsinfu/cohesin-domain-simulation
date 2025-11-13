#pragma once

#include <cstddef>
#include <iterator>


template<typename T>
class counting_iterator
{
public:
    using value_type = T;
    using reference = T;
    using pointer = void;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    explicit counting_iterator(value_type value)
    : _value(value)
    {
    }

    reference operator*() const
    {
        return _value;
    }

    counting_iterator& operator++()
    {
        ++_value;
        return *this;
    }

    counting_iterator operator++(int)
    {
        counting_iterator copy = *this;
        operator++();
        return copy;
    }

private:
    value_type _value;
};

template<typename T>
inline bool
operator==(counting_iterator<T> const& x, counting_iterator<T> const& y)
{
    return *x == *y;
}

template<typename T>
inline bool
operator!=(counting_iterator<T> const& x, counting_iterator<T> const& y)
{
    return !(x == y);
}
