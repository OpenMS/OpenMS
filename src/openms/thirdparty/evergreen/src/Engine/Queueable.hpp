#ifndef _QUEUEABLE_HPP
#define _QUEUEABLE_HPP

// Mixin to allow pointers of objects to be inserted into queues:
struct Queueable {
  // A member is faster than storing a map to which messages are in
  // queue in code:
  double priority;
  bool in_queue;

  Queueable():
    priority(0.0),
    in_queue(false)
  { }
};

#endif
