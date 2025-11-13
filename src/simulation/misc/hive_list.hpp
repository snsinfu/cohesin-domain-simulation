#pragma once

#include <cassert>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>


namespace detail
{
    template<typename From, typename To>
    using copy_const_t = std::conditional_t<std::is_const_v<From>, To const, To>;
}


// Unordered container with stable iterator after insertions and/or erasures.
template<typename T>
class hive_list
{
public:
    template<typename Container>
    class basic_iterator;

    using value_type      = T;
    using reference       = T&;
    using const_reference = T const&;
    using size_type       = std::size_t;
    using iterator        = basic_iterator<hive_list>;
    using const_iterator  = basic_iterator<hive_list const>;

    size_type size()  const noexcept
    {
        return _mask.size() - _free_slots.size();
    }

    bool empty() const noexcept
    {
        return size() == 0;
    }

    iterator       begin()        noexcept { return {this, 0}; }
    const_iterator begin()  const noexcept { return {this, 0}; }
    const_iterator cbegin() const noexcept { return {this, 0}; }
    iterator       end()        noexcept { return {this, _mask.size()}; }
    const_iterator end()  const noexcept { return {this, _mask.size()}; }
    const_iterator cend() const noexcept { return {this, _mask.size()}; }

    // Erases all elements in the container.
    void clear()
    {
        *this = {};
    }

    // Inserts an element with given value to an unspecified position in the
    // container. Returns an iterator pointing to the inserted element.
    iterator insert(value_type value)
    {
        auto const index = allocate_free_slot();
        if (index < _data.size()) {
            _data[index] = std::move(value);
            _mask[index] = 0;
        } else {
            _data.push_back(std::move(value));
            _mask.push_back(0);
        }
        return {this, index};
    }

    // Erases an element pointed to by the given iterator. Returns a value
    // moved out of the erased element.
    value_type erase(const_iterator const& it)
    {
        auto const index = it.cursor();
        assert(it.container() == this);
        assert(index < _mask.size());
        assert(_mask[index] == 0);
        value_type value = std::move(_data[index]);
        _mask[index] = 1;
        _free_slots.push_back(index);
        return value;
    }

    // Creates an iterator from an integral value rerutned by the
    // iterator::cursor() function.
    iterator locate(size_type cursor)
    {
        return {this, cursor};
    }

    const_iterator locate(size_type cursor) const
    {
        return {this, cursor};
    }

    // Creates an iterator from an object stored in this container.
    iterator locate_object(reference obj)
    {
        return locate(size_type(&obj - _data.data()));
    }

    iterator locate_object(const_reference obj)
    {
        return locate(size_type(&obj - _data.data()));
    }

    const_iterator locate_object(const_reference obj) const
    {
        return locate(size_type(&obj - _data.data()));
    }


private:
    size_type allocate_free_slot() noexcept
    {
        if (_free_slots.empty()) {
            return _mask.size();
        }
        auto const index = _free_slots.back();
        _free_slots.pop_back();
        assert(index < _mask.size());
        assert(_mask[index] == 1);
        return index;
    }

private:
    std::vector<T>         _data;
    std::vector<bool>      _mask;
    std::vector<size_type> _free_slots;
};


template<typename T>
template<typename Container>
class hive_list<T>::basic_iterator
{
public:
    using value_type        = T;
    using access_type       = detail::copy_const_t<Container, T>;
    using reference         = access_type&;
    using pointer           = access_type*;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;
    using size_type         = hive_list::size_type;

    basic_iterator() = default;

    basic_iterator(Container* hive, size_type index)
    : _hive{hive}, _index{index}
    {
        skip_masked_slots();
    }

    // Support iterator-to-const_iterator conversion.
    operator basic_iterator<Container const>() const noexcept
    {
        return {_hive, _index};
    }

    bool operator==(basic_iterator const& other) const
    {
        return _index == other._index;
    }

    bool operator!=(basic_iterator const& other) const
    {
        return !operator==(other);
    }

    reference operator*() const
    {
        return _hive->_data[_index];
    }

    pointer operator->() const
    {
        return &operator*();
    }

    basic_iterator& operator++()
    {
        _index++;
        skip_masked_slots();
        return *this;
    }

    basic_iterator operator++(int)
    {
        auto copy = *this;
        operator++();
        return copy;
    }

    // Returns the cursor value for the element pointed to by this iterator.
    //
    // A cursor is an integer associated with an element in a hive_list
    // container. Its value is no greater than the number of elements that have
    // been inserted to a container so far.
    //
    // The cursor is valid and may be passed to the locate() member function
    // of the associated hive_list instance until the element pointed to by
    // this iterator gets erased. So, a cursor can be used as a compact
    // alternative to an iterator.
    size_type cursor() const noexcept
    {
        return _index;
    }

    // Returns a pointer to the associated hive_list container object.
    Container* container() const noexcept
    {
        return _hive;
    }

private:
    void skip_masked_slots()
    {
        auto const& mask = _hive->_mask;
        while (_index < mask.size() && mask[_index]) {
            _index++;
        }
    }

private:
    Container* _hive  = nullptr;
    size_type  _index = 0;
};
