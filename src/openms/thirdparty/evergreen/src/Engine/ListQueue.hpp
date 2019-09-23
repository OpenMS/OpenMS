#ifndef _LISTQUEUE_HPP
#define _LISTQUEUE_HPP

#include <list>

template <typename VARIABLE_KEY>
class ListQueue {
protected:
  std::list<Edge<VARIABLE_KEY>* > _next_edges;
public:
  bool is_empty() const {
    return _next_edges.size() == 0;
  }
  std::size_t size() const {
    return _next_edges.size();
  }

  void push_if_not_in_queue(Edge<VARIABLE_KEY>* val) {
    if (val->in_queue)
      return;

    _next_edges.push_back(val);
    val->in_queue = true;
  }
  Edge<VARIABLE_KEY>* pop_next() {
    #ifdef ENGINE_CHECK
    assert( ! is_empty() );
    #endif

    Edge<VARIABLE_KEY>* res = _next_edges.front();
    _next_edges.pop_front();
    res->in_queue = false;
    return res;
  }
  void print(std::ostream & os) const {
    os << "Size " << size() << std::endl;

    for (const Edge<VARIABLE_KEY>* val : _next_edges) {
      os << val << " from " << val->source << " to " << val->dest << std::endl;
    }
    os << std::endl;
  }
};

#endif
