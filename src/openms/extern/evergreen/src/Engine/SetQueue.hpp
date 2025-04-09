#ifndef _SETQUEUE_HPP
#define _SETQUEUE_HPP

#include <set>
#include <unordered_map>
#include <iostream>

template <typename VARIABLE_KEY>
class SetQueue {
protected:
  double _max_priority;
  // This needs to use std::set instead of std::unordered_set because
  // the ordering is used (it used in a manner similar to a heap, so a
  // priority queue may also be an option).
  std::set<double> _priorities;

  // Uses unordered_set<T> since ordering is not used.
  std::unordered_map<double, std::unordered_set<Edge<VARIABLE_KEY>*> > _priorities_to_values;
  std::size_t _size;
public:
  SetQueue():
    _size(0)
  { }

  bool is_empty() const {
    return _size == 0;
  }
  std::size_t size() const {
    return _size;
  }

  double max_priority() const {
    #ifdef ENGINE_CHECK
    assert( ! is_empty() );
    #endif
    
    return _max_priority;
  }
  bool contains_priority(double priority) const {
    auto iter = _priorities.find(priority);
    return iter != _priorities.end();
  }

  void push_or_update(Edge<VARIABLE_KEY>* val, double new_priority) {
    if (val->in_queue)
      update_priority(val, new_priority);
    else {
      val->priority = new_priority;
      push(val);
    }
  }

  void push(Edge<VARIABLE_KEY>* val) {
    #ifdef ENGINE_CHECK
    assert(! val->in_queue);
    #endif

    if (! contains_priority(val->priority)) {
      _priorities.insert(val->priority);
      _priorities_to_values[val->priority] = std::unordered_set<Edge<VARIABLE_KEY>*>();
    }

    // Now both priorities and _priorities_to_values contain priority:
    std::unordered_set<Edge<VARIABLE_KEY>*> & vals_at_priority = _priorities_to_values[val->priority];
    #ifdef ENGINE_CHECK
    assert( vals_at_priority.find(val) == vals_at_priority.end() && "Value already in Queue");
    #endif
    vals_at_priority.insert(val);
    if (_size == 0 || val->priority > _max_priority)
      _max_priority = val->priority;
    ++_size;

    val->in_queue = true;
  }
  Edge<VARIABLE_KEY>* pop_max() {
    #ifdef ENGINE_CHECK
    assert( ! is_empty() );
    #endif

    double priority = max_priority();
    std::unordered_set<Edge<VARIABLE_KEY>*> & vals_at_priority = _priorities_to_values[priority];

    // Erase the max from the inner set:
    auto iter = vals_at_priority.begin();
    Edge<VARIABLE_KEY>* result = *iter;
    #ifdef ENGINE_CHECK
    assert(result->in_queue);
    #endif

    vals_at_priority.erase( iter );
    
    // If the inner set has become empty, remove the empty set and
    // remove the priority from the collection of all priorities:
    if (vals_at_priority.size() == 0) {
      _priorities_to_values.erase(priority);
      _priorities.erase(priority);
    }
    
    --_size;
    if ( ! is_empty() )
      _max_priority = *_priorities.rbegin();

    result->in_queue = false;
    return result;
  }
  void remove(Edge<VARIABLE_KEY>* val) {
    #ifdef ENGINE_CHECK
    assert(contains_priority(val->priority) && "Error: Priority to update not in queue");
    #endif
    // remove (possibly from middle):
    --_size;
    std::unordered_set<Edge<VARIABLE_KEY>*> & vals_at_priority = _priorities_to_values.find(val->priority)->second;
    #ifdef ENGINE_CHECK
    assert(vals_at_priority.count(val) && "Error: Value at requested priority not in queue");
    #endif
    vals_at_priority.erase( val );
    if ( vals_at_priority.size() == 0 ) {
      _priorities_to_values.erase(val->priority);
      _priorities.erase(val->priority);
    }

    if ( ! is_empty() )
      _max_priority = *_priorities.rbegin();

    val->in_queue = false;
  }
  void update_priority(Edge<VARIABLE_KEY>* val, double new_priority) {
    // Note: if new_priority < old_priority, a Fibonacci heap could do
    // this in amortized O(1):
    #ifdef ENGINE_CHECK
    assert(val->in_queue);
    #endif

    remove(val);
    val->priority = new_priority;
    push(val);
  }

  void print(std::ostream & os) const {
    os << "Size " << size() << std::endl;
    for (double priority : _priorities) {
      os << "Priority " << priority << " ";
      const std::unordered_set<Edge<VARIABLE_KEY>*> & vals_at_priority = _priorities_to_values.find(priority)->second;

      for (const Edge<VARIABLE_KEY>* val : vals_at_priority) {
	os << val;
	//	os << " " << val->get_possibly_outdated_message();
	os << " " << val->priority;
      }
      os << std::endl;
    }
  }
};

#endif
