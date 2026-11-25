// Index-based value iterator for element containers.
//
// nb::make_iterator stores raw C++ begin()/end() iterators, so mutating the
// container while a Python iteration is in flight dereferences invalidated
// memory (the push_back-reallocation failure class of issue #9792). This
// iterator holds a strong reference to the parent plus an index instead:
// every step bounds-checks against the container's *current* size and returns
// an owned copy of the element. Mutation during iteration is therefore
// defined the same way as for a Python list -- appended elements are visited,
// shrinking ends the loop early -- and can never touch freed memory.
//
// Requires the container to expose size() and operator[](index), or a
// specialization of IndexAccess for a container whose Python elements are
// counted or reached differently. Containers without index access (e.g.
// Param's tree iterator) should materialize a snapshot list instead.
#pragma once

#include <nanobind/nanobind.h>

namespace pyopenms_iter
{

// How the iterator counts and copies the elements of C
template <typename C>
struct IndexAccess
{
  static size_t size(C& c) { return c.size(); }
  static auto at(C& c, size_t i) { return c[i]; }  // by value: an owned copy
};

template <typename C>
struct IndexValueIterator
{
  nanobind::object owner;  // strong reference: keeps the container alive
  C* container = nullptr;  // stable address of the nanobind-owned instance
  size_t i = 0;
};

template <typename C>
IndexValueIterator<C> make_index_value_iterator(nanobind::object self)
{
  return IndexValueIterator<C>{self, &nanobind::cast<C&>(self), 0};
}

template <typename C>
void bind_index_value_iterator(nanobind::module_& m, const char* name)
{
  namespace nb = nanobind;
  using It = IndexValueIterator<C>;
  nb::class_<It>(m, name)
      .def("__iter__", [](nb::object self) { return self; })
      .def("__next__", [](It& it) {
        if (it.i >= IndexAccess<C>::size(*it.container)) { throw nb::stop_iteration(); }
        return IndexAccess<C>::at(*it.container, it.i++);  // by value: an owned copy per element
      });
}

}  // namespace pyopenms_iter
