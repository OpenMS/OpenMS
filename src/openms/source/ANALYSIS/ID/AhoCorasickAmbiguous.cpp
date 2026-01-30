// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <algorithm>
#include <cassert>
#include <queue>
#include <tuple>

namespace OpenMS
{
  /**
    @brief given an ambAA or the 'superAA' @p aaa, return a range of AA's which need to be explored.
  */
  inline constexpr std::tuple<AA,AA> _ambiguitiesOf(const AA aaa)
  {
    static_assert(++AA('D') == AA('N'));  // must be neighbors
    static_assert(++AA('E') == AA('Q'));  // must be neighbors
    static_assert(++AA('I') == AA('L'));  // must be neighbors
    // row of jump table must be B,J,Z,X,$:
    static_assert(++AA('B') == AA('J'));  // make sure the table is ordered as we expect
    static_assert(++AA('J') == AA('Z'));  // make sure the table is ordered as we expect
    static_assert(++AA('Z') == AA('X'));  // make sure the table is ordered as we expect
    static_assert(++AA('X') == AA('$'));  // make sure the table is ordered as we expect
    
    // jump table:                    start of scouts 
    //                                         end of scouts (including)
    constexpr const AA jump[5][2] = {{AA('D'), AA('N')},  // B = D,N
                                     {AA('I'), AA('L')},  // J = I,L
                                     {AA('E'), AA('Q')},  // Z = E,Q
                                     {AA('A'), AA('V')},  // X = A..V
                                     {AA('A'), AA('X')}}; // $ = A..X
    
    // which line of jump table do we need?
    const auto line = (aaa - AA('B'))();
    assert(aaa.isAmbiguous());
    return {jump[line][0], jump[line][1]};
  }

  ACTrie::ACTrie(uint32_t max_aaa, uint32_t max_mm) : max_aaa_(max_aaa), max_mm_(max_mm)
  { // create root node:
    trie_.emplace_back();
  }

  ACTrie::~ACTrie() = default;

  void ACTrie::addNeedle(const std::string& needle)
  {
    Index cn {0}; // start at root
    for (auto c : needle) // OMS_CODING_TEST_EXCLUDE
    {
      AA aa(c);
      // make sure invalid chars raise an exception
      if (aa.isValidForPeptide())
      {
        cn = add_(cn, aa);
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("Invalid amino acid"), std::string(1, c));
      }
    }
    // add hit to last node
    trie_[cn()].depth_and_hits.has_hit = 1;
    
    // remember a needle ends here
    if (vec_index2needles_.size() <= cn()) // make sure there is enough space
    { // increase to next power of 2
      vec_index2needles_.resize(std::bit_ceil(cn() + 1)); // +1 since bit_ceil(x)^2 == x iff x is a power of 2 (i.e. no resize takes place)
    }
    vec_index2needles_[cn()].push_back(needle_count_);
    ++needle_count_;
  }

  void ACTrie::addNeedles(const std::vector<std::string>& needles)
  {
    for (const auto& s : needles) 
    {
      addNeedle(s);
    }
  }

  void ACTrie::addNeedlesAndCompress(const std::vector<std::string>& needles)
  {
    for (const auto& s : needles)
    {
      addNeedle(s);
    }
    compressTrie();
  }

  void ACTrie::compressTrie()
  {
    // final BFS tree we want to create
    std::vector<ACNode> bfs_tree;
    bfs_tree.reserve(trie_.size());

    // translate old naive node index to new node index in BFS
    decltype(vec_index2needles_) bfs_index2_needles;
    bfs_index2_needles.resize(trie_.size());

    // points to the parent node for each node in the final BFS tree.
    // (needed for suffix construction)
    std::vector<Index> tmp_parents;
    tmp_parents.reserve(trie_.size());

    // contains nodes in breadth first search order
    std::queue<Index> bfs_q; 

    // lambda for each pop operation on the queue
    auto bfs_op = [&bfs_q, &bfs_tree, &bfs_index2_needles, &tmp_parents, this](Index current_index)
    {
      auto bfs_index = bfs_tree.size();
      // add current node to new trie
      bfs_tree.push_back(trie_[current_index()]);
      auto& bfs_node = bfs_tree.back();
      if (current_index() < vec_index2children_naive_.size() && vec_index2children_naive_[current_index()].size() != 0)
      {
        auto& children = vec_index2children_naive_[current_index()];
        bfs_node.nr_children = ACNode::ChildCountType(children.size());
        std::sort(children.begin(), children.end(), // sort children by edge label (so they have the same order as in the bitset)
                  [&](const Index& a, const Index& b) { return trie_[a()].edge() < trie_[b()].edge(); }); 

        for (const auto child : children)
        {
          bfs_node.children_bitset.set(trie_[child()].edge()); // set child exists
          bfs_q.push(child);
          tmp_parents.emplace_back(Index::T(bfs_index)); // the parent will be added at index = tmp_tree.size()
        }
      }

      bfs_index2_needles[bfs_index] = std::move(vec_index2needles_[current_index()]);
    };

    // create root manually
    tmp_parents.emplace_back(0); // root has no parents (points to itself)
    bfs_op(0);                // adds parents to 'tmp_parents' and children to 'BFS'
    ACNode& root = bfs_tree.back();
    root.first_child = 1; // we know the first child will start at index 1 (since root is index 0)
    while (!bfs_q.empty())
    {
      auto& last_node = bfs_tree.back(); // previous node in final tree
      Index current_index = bfs_q.front();
      bfs_q.pop();
      bfs_op(current_index);
      // update where to find the children of the current node
      // --> its where the children of the previous node end
      bfs_tree.back().first_child = last_node.first_child() + last_node.nr_children;
    }

    // switch to BFS trie
    trie_ = std::move(bfs_tree);
    vec_index2needles_.swap(bfs_index2_needles);

    /////////////////////////////////////////////////////////////
    // compute suffix links (could also be done while creating the trie, but it would make the code more complex)
    // .. and hit flag
    trie_[0].suffix = 0; // must point to itself

    /** 
    auto printTrie = [&]() {
      int prev_depth = 0;
      for (size_t i = 0; i < trie_.size(); ++i)
      {
        if (prev_depth != (int)trie_[i].depth_and_hits.depth)
        {
          prev_depth = (int)trie_[i].depth_and_hits.depth;
          std::cout << " ***\n";
        }
        std::cout << "node " << i << " (edge " << trie_[i].edge.toChar() << ", depth " << (int)trie_[i].depth_and_hits.depth << ", hit "
                  << (int)trie_[i].depth_and_hits.has_hit << ", suffix " << trie_[i].suffix() << ", first child " << (int)trie_[i].first_child()
                  << '\n';
      }
    };
    */

    /// first, build suffix links as usual (without path compression)
    // start at depth = 2, since setting suffix links and has_hit for depth 1 is not needed (lvl 1 already points to root)
    for (size_t i = 1 + (size_t)trie_[0].nr_children; i < trie_.size(); ++i)
    {
      Index parent = tmp_parents[i];
      trie_[i].suffix = follow_(trie_[parent()].suffix(), trie_[i].edge);
      trie_[i].depth_and_hits.has_hit |= trie_[trie_[i].suffix()].depth_and_hits.has_hit;
    }

    // second pass over trie, do path compression of suffix links where possible
    // start at depth = 2, since setting suffix links and has_hit for depth 1 is not needed (lvl 1 already points to root)
    size_t path_compression_count = 0;
    for (size_t i = 1 + (size_t)trie_[0].nr_children; i < trie_.size(); ++i)
    {
      const bool suffix_is_hit = trie_[trie_[i].suffix()].depth_and_hits.has_hit;
      if (suffix_is_hit) continue; // don't touch suffices which have hits (since we do not have separate output links)

      const auto children_node = trie_[i].children_bitset.bits;
      const auto childen_suffix = trie_[trie_[i].suffix()].children_bitset.bits;
      if ((children_node | childen_suffix) == children_node) // suffix' node has no extra children compared to this node `i`
      { 
        // point to suffix's suffix (path compression)
        trie_[i].suffix = trie_[trie_[i].suffix()].suffix;
        path_compression_count++;
      }
    }
    OPENMS_LOG_INFO << "ACTrie::compressTrie(): created BFS trie with " << trie_.size() << " nodes (path compression skipped " << path_compression_count
                    << " suffix nodes)\n";
    //printTrie();

    vec_index2children_naive_.clear(); // not needed anymore
  }

  size_t ACTrie::getNeedleCount() const
  {
    return needle_count_;
  }

  void ACTrie::setMaxAAACount(const uint32_t max_aaa)
  {
    max_aaa_ = max_aaa;
  }

  uint32_t ACTrie::getMaxAAACount() const
  {
    return max_aaa_;
  }

  void ACTrie::setMaxMMCount(const uint32_t max_mm)
  {
    max_mm_ = max_mm;
  }

  uint32_t ACTrie::getMaxMMCount() const
  {
    return max_mm_;
  }

  bool ACTrie::nextHits(ACTrieState& state) const
  {
    state.hits.clear();
    assert(vec_index2children_naive_.empty()); // make sure compressTrie was called
    nextHitsNoClear_(state);
    return !state.hits.empty();
  }

  void ACTrie::getAllHits(ACTrieState& state) const
  {
    state.hits.clear();
    assert(vec_index2children_naive_.empty()); // make sure compressTrie was called
    while (nextHitsNoClear_(state)) {};
  }

  bool ACTrie::nextHitsNoClear_(ACTrieState& state) const
  {
    std::vector<Hit>& hits = state.hits;
    for (AA aa = state.nextValidAA(); aa.isValid(); aa = state.nextValidAA())
    {
      state.tree_pos = stepPrimary_(state.tree_pos, aa, state);
      // deal with scouts in queue: doing it now (instead of after the primary ends) benefits from hot caches
      // and a lot less memory (since only hits from current scouts are found)
      while (!state.scouts.empty())
      {
        ACScout& sp = state.scouts.front();
        // let scout traverse the tree until it dies. This might add new scouts to the queue.
        while (stepScout_(sp, state));
        state.scouts.pop();
      }
      if (addHits_(state.tree_pos, state.textPos(), hits))
      {
        return true;
      };
    }



    return false;
  }

  Index ACTrie::add_(const Index index, const AA label)
  {
    if (vec_index2children_naive_.size() <= index())
    { // doubling is not enough
      vec_index2children_naive_.resize(std::bit_ceil(index() + 1));  // +1 since bit_ceil(x)^2 == x iff x is a power of 2 (i.e. no resize takes place)
    }
    Index ch = findChildNaive_(index, label);
    if (ch.isInvalid())
    {
      // remember index of new node we are about to create
      ch.pos() = Index::T(trie_.size());
      // create new node with label and depth
      trie_.emplace_back(label, trie_[index()].depth_and_hits.depth + 1);
      // add child to parent

      vec_index2children_naive_[index()].push_back(ch);
    }
    return ch;
  }


  bool ACTrie::addHits_(Index i, const size_t text_pos, std::vector<Hit>& hits) const
  {
    size_t hits_before = hits.size();
    // hits from current node; return true if going upstream has more hits..
    auto collect = [&]() {
      if (trie_[i()].depth_and_hits.has_hit)
      {
        const auto needle_length = trie_[i()].depth_and_hits.depth;
        const auto text_start = text_pos - needle_length;
        for (const auto needle_idx : vec_index2needles_.at(i()))
        {
          hits.emplace_back(needle_idx, needle_length, Hit::T(text_start));
        }
        return true;
      }
      return false;
    };

    // follow chain of suffix nodes until a node does not have hits anymore
    while (collect())
    {
      i = trie_[i()].suffix;
    }
    return hits_before != hits.size();
  }

  bool ACTrie::addHitsScout_(Index i, const ACScout& scout, const size_t text_pos, std::vector<Hit>& hits, const int current_scout_depths) const
  {
    size_t hits_before = hits.size();
    // hits from current node; return true if going upstream has more hits..
    auto collect = [&]() {
      if (trie_[i()].depth_and_hits.has_hit)
      {
        const auto needle_length = trie_[i()].depth_and_hits.depth;
        const auto text_start = text_pos - needle_length;
        // we want the first AAA of the scout to be part of the hit; otherwise that hit will be reported by shorter sub-scouts or the Primary
        if (current_scout_depths - needle_length >= scout.max_prefix_loss_leftover) 
        {
          return false;
        }
        for (const auto needle_idx : vec_index2needles_.at(i()))
        {
          hits.emplace_back(needle_idx, needle_length, Hit::T(text_start));
        }
        return true;
      }
      return false;
    };

    // follow chain of suffix nodes until a node does not have hits anymore
    while (collect())
    {
      i = trie_[i()].suffix;
    }
    return hits_before != hits.size();
  }
  
  Index ACTrie::follow_(const Index i, const AA aa) const
  {
    Index ch = findChildBFS_(i, aa);
    // has direct child (could also be an ambiguous AA - we don't care as long as a needle did contain that character)
    if (ch.isValid())
    {
      return ch;
    }

    // no direct child; are we at root?
    if (i.pos() == 0)
    {
      return 0;
    }
      
    // follow from the suffix...
    Index suf = trie_[i.pos()].suffix;
    assert(suf.isValid());
    return follow_(suf, aa);
  }

  bool ACTrie::followScout_(ACScout& scout, const AA edge, ACTrieState& state) const
  {
    // let scout follow the original edge
    Index j = follow_(scout.tree_pos, edge);
    const int new_depth = int(trie_[j()].depth_and_hits.depth);
    // did we loose a prefix?              old-depth                         new depth
    const int up_count = int(trie_[scout.tree_pos()].depth_and_hits.depth) - new_depth + 1;
    if (up_count >= scout.max_prefix_loss_leftover)
    { // scout is dead because it lost its AAA/MM
      return false;
    }
    // update the prefix length
    scout.max_prefix_loss_leftover -= up_count;
    scout.tree_pos = j;
    addHitsScout_(j, scout, scout.textPos(state), state.hits, new_depth);
    return true;
  }

  Index ACTrie::stepPrimary_(const Index i, const AA edge, ACTrieState& state) const
  {
    const bool consider_ambAA = max_aaa_ != 0;
    const bool consider_MM = max_mm_ != 0;
    
    // AAA
    if (edge.isAmbiguous())
    { // create AAA scouts?
      AA from(edge), to(edge);
      if (consider_ambAA)
      { // first try AAA's (since edge is AAA)
        std::tie(from, to) = _ambiguitiesOf(edge); // e.g. [D,N] for B; i.e. create two scouts; Primary will follow 'B' (if exists in pattern; if not, it will end up in root)
        createScouts_(i, from, to, state, max_aaa_ - 1, max_mm_);
      }
      // test all other AA's for mismatch
      if (consider_MM)
      { 
        createMMScouts_(i, from, to, edge, state, max_aaa_, max_mm_ - 1); // try a MM for all AA's other than [from...to]
      }
    }
    // edge is unambiguous:
    else if (consider_MM)
    { // try a MM for all AA's other than 'edge'
      createMMScouts_(i, edge, edge, edge, state, max_aaa_, max_mm_ - 1);
    }
  
    // Primary continues with the AA, no matter what it was...
    Index ch = findChildBFS_(i, edge);

    // has direct child (could also be an ambiguous AA - we don't care as long as a needle did contain that character)
    if (ch.isValid())
    {
      return ch;
    }
    // are we at root?
    if (i() == 0)
    {
      return i;
    }
    // follow from the suffix...
    Index suf = trie_[i()].suffix;
    assert(suf.isValid());
    return follow_(suf, edge);
  }

  bool ACTrie::stepScout_(ACScout& scout, ACTrieState& state) const
  {
    for (AA edge = scout.nextValidAA(); edge.isValid(); edge = scout.nextValidAA())
    {
      const bool consider_ambAA = scout.max_aaa_leftover > 0;
      const bool consider_MM = scout.max_mm_leftover > 0;

      // AAA
      if (edge.isAmbiguous())
      { // create scouts from this scout?
        AA from(edge), to(edge);
        if (consider_ambAA)
        { // first try AAA's (since edge is AAA)
          std::tie(from, to) = _ambiguitiesOf(edge);
          ACScout sp_temp = scout;
          --sp_temp.max_aaa_leftover;
          createSubScouts_(sp_temp, from, to, state);
        }
        // test all other superAA's for mismatch (except for AAA range, and the original edge itself)
        if (consider_MM)
        {
          ACScout sp_temp = scout;
          --sp_temp.max_mm_leftover;
          createMMSubScouts_(sp_temp, from, to, edge, state);
        }
      }
      else if (consider_MM) // edge is unambiguous
      { // try a MM for all superAA's other than 'edge'
        ACScout sp_temp = scout;
        --sp_temp.max_mm_leftover;
        createMMSubScouts_(sp_temp, edge, edge, edge, state);
      }

      // process the scout itself
      if (!followScout_(scout, edge, state)) return false;
    }
    return false; // end of query reached
  }

  void ACTrie::createMMScouts_(const Index i, const AA except_fromAA, const AA except_toAA, const AA except_edge, ACTrieState& state, const uint32_t aaa_left, const uint32_t mm_left) const
  {
    // create super-AA range, i.e. including the ambiguous AA's, since a peptide could contain an 'X', which we would like to match
    auto [from, to] = _ambiguitiesOf(AA('$'));
    for (AA mm_aa = from; mm_aa <= to; ++mm_aa)
    {
      if (mm_aa == except_fromAA)
      { // ignore this range
        mm_aa = except_toAA;
        continue;
      }
      // ignore edge from scout
      if (mm_aa == except_edge) 
      { 
        continue;
      }
      createScouts_(i, mm_aa, mm_aa, state, aaa_left, mm_left);
    }
  }

  void ACTrie::createMMSubScouts_(const ACScout& prototype, const AA except_fromAA, const AA except_toAA, const AA except_edge, ACTrieState& state) const
  {
    // create super-AA range, i.e. including the ambiguous AA's, since a peptide could contain an 'X', which we would like to match
    auto [from, to] = _ambiguitiesOf(AA('$'));
    for (AA mm_aa = from; mm_aa <= to; ++mm_aa)
    {
      if (mm_aa == except_fromAA)
      { // ignore this range
        mm_aa = except_toAA;
        continue;
      }
      // ignore edge from scout
      if (mm_aa == except_edge)
      {
        continue;
      }
      createSubScouts_(prototype, mm_aa, mm_aa, state);
    }
  }

  void ACTrie::createScouts_(const Index i, const AA fromAA, const AA toAA, ACTrieState& state, const uint32_t aaa_left, const uint32_t mm_left) const
  {
    for (AA aa = fromAA; aa <= toAA; ++aa)
    {
      Index scout_pos = follow_(i, aa); // call this using naive follow_(), which matches the exact char
      if (scout_pos() > 0) // not at root
      {
        const uint8_t depth = trie_[scout_pos()].depth_and_hits.depth;
        auto new_scout = state.scouts.emplace(state.textPosIt(), // the master already points to the next AA, so scout can start there
                                               scout_pos,
                                               aaa_left,
                                               mm_left,
                                               depth);
        // we might have found a hit already: report it
        addHits_(scout_pos, new_scout.textPos(state), state.hits);
      }
    }
  }

  void ACTrie::createSubScouts_(const ACScout& prototype, const AA fromAA, const AA toAA, ACTrieState& state) const
  {
    for (AA aa = fromAA; aa <= toAA; ++aa)
    {
      ACScout s(prototype);
      if (followScout_(s, aa, state))
      { // scout survived following the edge
        state.scouts.push(std::move(s));
      }
    }
  }

  Index ACTrie::findChildNaive_(Index parent, AA child_label)
  {
    for (Index child : vec_index2children_naive_[parent()])  // only a 4byte type: copy it
    {
      if (trie_[child.pos()].edge == child_label)
        return child;
    }
    return Index {};
  }

  /// Count how many bits are set in the bitset up to (not including) position i
  int countSetBitsUpTo(Bitset bs, unsigned int i)
  {
    static_assert(sizeof(bs) == 4, "Bitset must be 32 bits wide");
    if (0 == i)
    {
      return 0; // no bits to count
    }
    unsigned int high_bits_to_eliminate = 32 - i; // i should be <= 31 (we don't check for performance reasons)
    bs <<= (high_bits_to_eliminate);
    return bs.pop_count();
  }

  Index ACTrie::findChildBFS_(const Index parent, const AA child_label) const
  {
    // check if it exists
    if (trie_[parent()].children_bitset.test(child_label()) == 0)
    {
      return Index{}; // return invalid index
    }
    // child exists:
    // .. find its offset (children are ordered by edge label)
    auto child_offset = countSetBitsUpTo(trie_[parent()].children_bitset, child_label()); 
    return trie_[parent()].first_child() + child_offset;
  }

  void ACTrieState::setQuery(const std::string& haystack)
  {
    hits.clear();
    query_ = haystack;
    it_q_ = &query_[0];
    tree_pos = 0;
    while (!scouts.empty())
    {
      scouts.pop();
    }
  }

  size_t ACTrieState::textPos() const
  {
    return std::distance(&query_[0], it_q_);
  }

  const char* ACTrieState::textPosIt() const
  {
    return it_q_;
  }

  /// The current query

  const std::string& ACTrieState::getQuery() const
  {
    return query_;
  }

  AA ACTrieState::nextValidAA()
  {
    return OpenMS::nextValidAA(it_q_);
  }

  ACScout::ACScout(const char* query_pos, Index tree_pos, uint8_t max_aa, uint8_t max_mm, uint8_t max_prefix_loss) :
      it_query(query_pos), tree_pos(tree_pos), max_aaa_leftover(max_aa), max_mm_leftover(max_mm), max_prefix_loss_leftover(max_prefix_loss)
  {
  }

  size_t ACScout::textPos(const ACTrieState& state) const
  {
    return std::distance(&state.query_[0], it_query);
  }

  AA ACScout::nextValidAA()
  {
    return OpenMS::nextValidAA(it_query);
  }

  AA nextValidAA(const char*& it_q)
  {
    const char* it_q_local = it_q; // local copy; huge performance loss to work on it_q directly (due to reference)
    AA res {'?'}; // invalid
    while (*it_q_local != '\0')
    {
      res = AA(*it_q_local);
      ++it_q_local;
      if (res.isValid())
      {
        break;
      }
    }
    it_q = it_q_local; // update the reference to iterator
    return res;
  }

} // namespace OpenMS
