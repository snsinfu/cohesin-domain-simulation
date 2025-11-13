#pragma once

#include <cstddef>


template<typename Iterator>
class sized_iterator_range
{
public:
    using iterator  = Iterator;
    using size_type = std::size_t;
 
    inline sized_iterator_range(iterator begin, iterator end, size_type size)
    : _begin{begin}, _end{end}, _size{size}
    {
    }

    inline size_type size() const
    {
        return _size;
    }

    inline iterator begin() const
    {
        return _begin;
    }

    inline iterator end() const
    {
        return _end;
    }

private:
    iterator  _begin, _end;
    size_type _size;
};
