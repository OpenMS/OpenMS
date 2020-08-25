#ifndef _CONVOLUTIONTREE_HPP
#define _CONVOLUTIONTREE_HPP

#include "../PMF/PMF.hpp"
#include <limits>

// A convolution tree optimized for online processing. Note that the
// messages received through a given channel should never grow in
// support (otherwise, the cache of possible supports can become
// corrupted). This should always be the case in loopy belief
// propagation, since a growing product of PMFs will be passed. This
// is currently not checked because only the narrowed prior and
// narrowed likelihood are stored (not the raw prior and likelihood
// that were sent in). It could be checked, but at the expense of
// storing additional PMFs.

// Note: the recursive design of the convolution tree isn't the
// fastest or most robust (iterative would be better), but it easily
// allows lazy updating (i.e., messages are not propagated until a
// message out is requested). It could also be beneficial to allocate
// nodes in a single block rather, so that there is greater
// localization.

// The #define DISABLE_TRIM turns off trimming; since trimming has no
// real disadvantages, this is only used to measure the benefits of
// trimming.

class TreeNode {
#ifdef CONVOLUTIONTREE_CONVOLUTION_SIZE_CHECK
public:
  // For testing that supports are properly trimmed:
  static unsigned long largest_convolution_size;
#endif

protected:
  PMF _prior;
  PMF _likelihood;

  // Used to conservatively bound the support (this allows PMFs to be
  // narrowed in certain cases, which can dramatically improve
  // runtime). Note that the updating scheme updates these supports
  // before updating the related distribution; this is necessary to
  // allow the distribution to be fully trimmed (consider a tree where
  // nothing is yet cached, but where all inputs are binary and all
  // outputs are binary; without first setting support for all nodes
  // to binary, requesting the prior message out from the root will
  // produce long convolutions throughout the tree). This does have a
  // quirk where the supports will not be "dirtied" when a
  // distribution is narrowed additionally because initial narrowing
  // produces a distribution whose initial bounding box is full of
  // zeros. However, this case is addressed in the direction messages
  // have been requested, because changing the distrubution also
  // narrows both the support and distribution to the intersection of
  // both. The only case not addressed is dirtying the supports going
  // in the opposite direction. As a result, messages requested in
  // that direction may not benefit from this extra narrowing of
  // supports; however, dirtying them would be excessive, since this
  // case of a zero bounding box can be considered non-general, and
  // aggressively dirtying the supports in every direction could harm
  // the general runtime by making updating of supports no longer lazy
  // (as implemented, the dirtying cost can be amortized out).
  Vector<long> _minimum_possible_first_support;
  Vector<long> _maximum_possible_last_support;
  
  bool _prior_ready;
  bool _likelihood_ready;

  bool _support_from_below_ready;
  bool _support_from_above_ready;

  // Note that (with the exception of the root) sibling should exist,
  // and this function should never return NULL, because the tree
  // should be full.
  TreeNode* sibling_ptr() {
    if (parent->child_lhs == this)
      return parent->child_rhs;
    return parent->child_lhs;
  }
  bool has_children() {
    // Note that nodes have either 0 or 2 children since it is a full
    // binary tree; therefore, the following line can be simplified:
    //return child_lhs != NULL || child_rhs != NULL;
    return child_lhs != NULL;
  }
  void set_dependents_up_not_ready() {
    // Continue until reaching node that is not ready from below
    // (_prior_ready and _support_from_below_ready should match, but
    // use both to be tidy):
    if ( _prior_ready || _support_from_below_ready ) {
      _prior_ready = false;
      _support_from_below_ready = false;
      
      if (parent != NULL) {
	parent->set_dependents_up_not_ready();
	
	// If parent != NULL, sib should exist because the tree should
	// be full.
	TreeNode*sib = sibling_ptr();
	sib->set_dependents_down_not_ready();
      }
    }
  }
  void set_dependents_down_not_ready() {
    // Continue until reaching node that is not ready from above
    // (_likelihood_ready and _support_from_above_ready should match,
    // but use both to be tidy):
    if ( _likelihood_ready || _support_from_above_ready ) {
      _likelihood_ready = false;
      _support_from_above_ready = false;
      
      if (child_lhs != NULL)
	child_lhs->set_dependents_down_not_ready();
      if (child_rhs != NULL)
	child_rhs->set_dependents_down_not_ready();
    }
  }
  // add and sub return only the initialized sum or difference,
  // otherwise only add/subtract the initialized argument:
  inline static PMF add(const PMF & lhs, const PMF & rhs, double p) {
    if (lhs.dimension() == 0)
      return rhs;
    if (rhs.dimension() == 0)
      return lhs;

    #ifdef CONVOLUTIONTREE_CONVOLUTION_SIZE_CHECK  
    // Update the size of the maximum convolution performed:
    unsigned long n = std::max(lhs.table().flat_size(), rhs.table().flat_size());
    largest_convolution_size = std::max(largest_convolution_size, n);
    #endif

    return p_add(lhs, rhs, p);
  }
  inline static PMF sub(const PMF & lhs, const PMF & rhs, double p) {
    if (lhs.dimension() == 0)
      return rhs;
    if (rhs.dimension() == 0)
      return lhs;

    #ifdef CONVOLUTIONTREE_CONVOLUTION_SIZE_CHECK  
    // Update the size of the maximum convolution performed:
    unsigned long n = std::max(lhs.table().flat_size(), rhs.table().flat_size());
    largest_convolution_size = std::max(largest_convolution_size, n);
    #endif

    return p_sub(lhs, rhs, p);
  }
  void narrow_support_with(PMF & dist) {
    if (dist.dimension() != 0) {
      // Narrow dist to minimum and maximum supports:
      dist.narrow_support(_minimum_possible_first_support, _maximum_possible_last_support);
      
      // Narrow minimum and maximum supports to dist:
      for (unsigned char i=0; i<_minimum_possible_first_support.size(); ++i) {
	_minimum_possible_first_support[i] = std::max(_minimum_possible_first_support[i], dist.first_support()[i]);
	_maximum_possible_last_support[i] = std::min(_maximum_possible_last_support[i], long(dist.first_support()[i] + dist.table().view_shape()[i]) - 1);
      }
    }
  }
  void narrow_all() {
    #ifndef DISABLE_TRIM
    narrow_support_with(_likelihood);
    narrow_support_with(_prior);
    // Just in case prior narrows min/max supports, propagate that
    // change to likelihood:
    narrow_support_with(_likelihood);
    #endif
  }
  void update_prior(double p) {
    
    if ( ! _prior_ready ) {
      // Full binary tree means both must be non-null
      // simultaneously, so no need to consider cases where exactly
      // one is NULL:
      if (has_children()) {
	child_lhs->update_prior(p);
	child_rhs->update_prior(p);

	if (child_lhs->_prior_ready && child_rhs->_prior_ready)
	  set_prior( add(child_lhs->get_prior(p), child_rhs->get_prior(p), p) );
      }
    }
  }
  void update_likelihood(double p) {
    
    if ( ! _likelihood_ready ) {
      if (parent != NULL) {
	parent->update_likelihood(p);
	TreeNode*sib = sibling_ptr();
	sib->update_prior(p);
	if (parent->_likelihood_ready && sib->_prior_ready)
	  set_likelihood( sub(parent->get_likelihood(p), sib->get_prior(p), p) );
      }
    }
  }
  void update_support_from_below() {
    
    if ( ! _support_from_below_ready ) {
      if (child_lhs != NULL && child_rhs != NULL) {
	child_lhs->update_support_from_below();
	child_rhs->update_support_from_below();

	if ( child_lhs->_support_from_below_ready && child_rhs->_support_from_below_ready ) {
	  for (unsigned char i=0; i<_minimum_possible_first_support.size(); ++i) {
	    _minimum_possible_first_support[i] = std::max(_minimum_possible_first_support[i], child_lhs->_minimum_possible_first_support[i] + child_rhs->_minimum_possible_first_support[i]);
	    _maximum_possible_last_support[i] = std::min(_maximum_possible_last_support[i], child_lhs->_maximum_possible_last_support[i] + child_rhs->_maximum_possible_last_support[i]);
	  }

	  narrow_all();
	  _support_from_below_ready = true;
	}
      }
    }
  }

  void update_support_from_above() {
    
    if ( ! _support_from_above_ready ) {
      if (parent != NULL) {
	parent->update_support_from_above();
	TreeNode*sib = sibling_ptr();
	sib->update_support_from_below();

	if (parent->_support_from_above_ready && sib->_support_from_below_ready) {
	  // Note: This can be done more efficiently by inlining the
	  // following two lines into the loop below (memory
	  // allocation is more expensive than performing
	  // elementwise).
	  Vector<long> likelihood_minimum_possible_first_support = parent->_minimum_possible_first_support - sib->_maximum_possible_last_support;
	  Vector<long> likelihood_maximum_possible_last_support = parent->_maximum_possible_last_support - sib->_minimum_possible_first_support;

	  for (unsigned char i=0; i<likelihood_minimum_possible_first_support.size(); ++i) {
	    _minimum_possible_first_support[i] = std::max(_minimum_possible_first_support[i], likelihood_minimum_possible_first_support[i]);
	    _maximum_possible_last_support[i] = std::min(_maximum_possible_last_support[i], likelihood_maximum_possible_last_support[i]);
	  }

	  narrow_all();
	  _support_from_above_ready = true;
	}
      }
    }
  }
public:
  TreeNode *parent, *child_lhs, *child_rhs;
  
  TreeNode(unsigned char dimension):
    _minimum_possible_first_support(dimension),
    _maximum_possible_last_support(dimension),
    _prior_ready(false),
    _likelihood_ready(false),
    _support_from_below_ready(false),
    _support_from_above_ready(false),
    parent(NULL),
    child_lhs(NULL),
    child_rhs(NULL)
  {
    for (unsigned char i=0; i<dimension; ++i) {
      _minimum_possible_first_support[i] = std::numeric_limits<long>::min();
      _maximum_possible_last_support[i] = std::numeric_limits<long>::max();
    }
  }
  void set_prior(PMF && pmf) {
    _prior = std::move(pmf);
    narrow_all();

    // Set all dependents of prior to be unready:

    // Note: this will be called multiple times when updating many
    // nodes that are on a path to the root, but it will not matter
    // since the runtime will be O(1) per call after the first call:
    set_dependents_up_not_ready();

    // Set local prior to be ready:
    _prior_ready = true;

    if (! has_children())
      _support_from_below_ready = true;
  }
  void set_likelihood(PMF && pmf) {
    
    _likelihood = std::move(pmf);
    narrow_all();

    // Set all dependents of likelihood to be unready:

    // Note: this will be called multiple times when updating many
    // nodes that are on a path to the root, but it will not matter
    // since the runtime will be O(1) per call after the first call:
    set_dependents_down_not_ready();

    // Set local likelihood to be ready:
    _likelihood_ready = true;

    if (parent == NULL)
      _support_from_above_ready = true;
  }
  const PMF & get_prior(double p) {
    update_support_from_above();
    update_prior(p);
    #ifdef ENGINE_CHECK
    assert(_prior_ready);
    #endif
    return _prior;
  }
  const PMF & get_likelihood(double p) {
    update_support_from_above();
    update_likelihood(p);
    #ifdef ENGINE_CHECK
    assert(_likelihood_ready);
    #endif
    return _likelihood;
  }
  void add_child_lhs(TreeNode*lhs) {
    child_lhs = lhs;
    lhs->parent = this;
  }
  void add_child_rhs(TreeNode*rhs) {
    child_rhs = rhs;
    rhs->parent = this;
  }
  // For debugging:
  void print(std::ostream & os, unsigned int depth=0) {
    for (unsigned int i=0; i<3*depth; ++i)
      os << " ";
    os << this << " prior&support " << _prior_ready << _support_from_below_ready << " likelihood&support " << _likelihood_ready << _support_from_above_ready << " min/max possible support " << _minimum_possible_first_support << " " << _maximum_possible_last_support << " prior/likelihood " << _prior << " " << _likelihood << std::endl;
    
    if (child_lhs != NULL && child_rhs != NULL) {
      child_lhs->print(os, depth+1);
      child_rhs->print(os, depth+1);
    }
  }
};

#ifdef CONVOLUTIONTREE_CONVOLUTION_SIZE_CHECK  
unsigned long TreeNode::largest_convolution_size = 0;
#endif

class ConvolutionTree {
protected:
  const unsigned char _dimension;
  const double _p;

  TreeNode * _root;
  std::vector<TreeNode*> _inputs;
  // _output is the same as _root:

  // Construct a full binary tree with n leaves:
  TreeNode* create_tree(unsigned long number_priors_to_add) {
    TreeNode*res = new TreeNode(_dimension);
    if (number_priors_to_add > 1) {
      // When n == 1, allocate a single leaf. Otherwise, allocate
      // floor(n/2) leaves on the left subtree and n - floor(n/2)
      // leaves on the right subtree. When n>1, the left subtree will
      // get at least one child and the right subtree will get at
      // least one child, guaranteeing a full tree (when n != 1,
      // two subtrees will be created).
      res->add_child_lhs( create_tree(number_priors_to_add >> 1) );
      res->add_child_rhs( create_tree( number_priors_to_add - (number_priors_to_add >> 1) ) );
    }
    else {
      // leaf: add it to _inputs:
      _inputs.push_back(res);
    }
    return res;
  }

  void destroy_tree(TreeNode*&node) {
    if (node == NULL)
      return;

    if (node->child_lhs != NULL)
      destroy_tree(node->child_lhs);
    if (node->child_rhs != NULL)
      destroy_tree(node->child_rhs);
    delete node;
    node = NULL;
  }
  
public:

  // For debugging:
  void print(std::ostream & os) {
    _root->print(os);
  }

  ConvolutionTree(unsigned long number_priors_to_add, const unsigned char dim, const double p):
    _dimension(dim),
    _p(p)
  {
    _root = create_tree(number_priors_to_add);
  }

  ~ConvolutionTree() {
    destroy_tree(_root);
  }

  void receive_message_in(unsigned long index, PMF msg) {
    if (index < _inputs.size())
      // If the index is in 0, ... n-1, it refers to an input prior:
      _inputs[index]->set_prior(std::move(msg));
    else
      // Otherwise, the index refers to the output likelihood:
      _root->set_likelihood(std::move(msg));
  }

  PMF get_message_out(unsigned long index) {
    // The check as to whether this message can be computed will be
    // handled by the TreeNode types.
    
    if (index < _inputs.size()) {
      // If the index is in 0, ... n-1, it refers to an input prior:
      return _inputs[index]->get_likelihood(_p);
    }
    
    // Otherwise, the index refers to the output likelihood:
    return _root->get_prior(_p);
  }

  unsigned char dimension() const {
    return _dimension;
  }
};

#endif
