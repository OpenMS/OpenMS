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

#include <cassert>
#include <queue>

namespace OpenMS
{
  /**
    @brief given an ambAA or the 'superAA' @p aa, return a range of AA's which need to be spawned.
  */
  inline constexpr std::tuple<AA,AA> _getSpawnRange(const AA aa)
  {
    static_assert(++AA('D') == AA('N'));  // must be neighbors
    static_assert(++AA('E') == AA('Q'));  // must be neighbors
    static_assert(++AA('I') == AA('L'));  // must be neighbors
    // row of jump table must be B,J,Z,X,$:
    static_assert(++AA('B') == AA('J'));  // make sure the table is ordered as we expect
    static_assert(++AA('J') == AA('Z'));  // make sure the table is ordered as we expect
    static_assert(++AA('Z') == AA('X'));  // make sure the table is ordered as we expect
    static_assert(++AA('X') == AA('$'));  // make sure the table is ordered as we expect
    
    // jump table:                    start of spawns 
    //                                         end of spawns (including)
    constexpr const AA jump[5][2] = {{AA('D'), AA('N')},  // B = D,N
                                     {AA('I'), AA('L')},  // J = I,L
                                     {AA('E'), AA('Q')},  // Z = E,Q
                                     {AA('A'), AA('V')},  // X = A..V
                                     {AA('A'), AA('X')}}; // $ = A..X
    
    // which line of jump table do we need?
    const auto line = (aa - AA('B'))();
    assert(aa.isAmbiguous());
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
    umap_index2needles_[cn()].push_back(needle_count_);
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
    decltype(umap_index2needles_) bfs_index2_needles;

    // points to the parent node for each node in the final BFS tree.
    // (needed for suffix construction)
    std::vector<Index> tmp_parents;
    tmp_parents.reserve(trie_.size());

    // contains nodes in breadth first search order
    std::queue<Index> bfs_q; 

    // lambda for each pop operation on the queue
    auto bfs_op = [&bfs_q, &bfs_tree, &bfs_index2_needles, &tmp_parents, this](Index current_index) {
      // add children to BFS
      const auto& children = umap_index2children_naive_[current_index];
      auto bfs_index = bfs_tree.size();
      for (const auto& child : children)
      {
        bfs_q.push(child);
        tmp_parents.emplace_back(Index::T(bfs_index)); // the parent will be added at index = tmp_tree.size()
      }
      // add current node to new trie
      bfs_tree.push_back(trie_[current_index()]);
      bfs_tree.back().nr_children = ACNode::ChildCountType(children.size());

      bfs_index2_needles[Index::T(bfs_index)] = std::move(umap_index2needles_[current_index()]);
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
    umap_index2needles_ = std::move(bfs_index2_needles);

    // compute suffix links (could also be done while creating the trie, but it would make the code more complex)
    // .. and hit flag
    trie_[0].suffix = 0; // must point to itself

    Index old_parent {0};
    // start at depth = 2, since setting suffix links and has_hit for depth 1 is not needed
    // and will lead to suffix links pointing to itself
    for (size_t i = 1 + (size_t)trie_[0].nr_children; i < trie_.size(); ++i)
    {
      Index parent = tmp_parents[i];
      trie_[i].suffix = follow_(trie_[parent()].suffix(), trie_[i].edge);
      trie_[i].depth_and_hits.has_hit |= trie_[trie_[i].suffix()].depth_and_hits.has_hit;
    }
    umap_index2children_naive_.clear(); // not needed anymore
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
    assert(umap_index2children_naive_.empty()); // make sure compressTrie was called
    nextHitsNoClear_(state);
    return !state.hits.empty();
  }

  void ACTrie::getAllHits(ACTrieState& state) const
  {
    state.hits.clear();
    assert(umap_index2children_naive_.empty()); // make sure compressTrie was called
    while (nextHitsNoClear_(state)) {};
  }

  bool ACTrie::nextHitsNoClear_(ACTrieState& state) const
  {
    std::vector<Hit>& hits = state.hits;
    for (AA aa = state.nextValidAA(); aa.isValid(); aa = state.nextValidAA())
    {
      state.tree_pos = stepMaster_(state.tree_pos, aa, state);
      if (addHits_(state.tree_pos, state.textPos(), hits))
      {
        return true;
      };
    }

    // deal with spawns in queue
    while (!state.spawns.empty())
    {
      ACSpawn& sp = state.spawns.front();
      // let spawn traverse the tree until it dies. This might add new spawns to the queue.
      while (stepSpawn_(sp, state));
      state.spawns.pop();
    }

    return false;
  }

  Index ACTrie::add_(const Index index, const AA label)
  {
    Index ch = findChildNaive_(index, label);
    if (ch.isInvalid())
    {
      // remember index of new node we are about to create
      ch.pos() = Index::T(trie_.size());
      // create new node with label and depth
      trie_.emplace_back(label, trie_[index()].depth_and_hits.depth + 1);
      // add child to parent
      umap_index2children_naive_[index].push_back(ch);
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
        for (const auto needle_idx : umap_index2needles_.at(i()))
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

  bool ACTrie::addHitsSpawn_(Index i, const ACSpawn& spawn, const size_t text_pos, std::vector<Hit>& hits, const int current_spawn_depths) const
  {
    size_t hits_before = hits.size();
    // hits from current node; return true if going upstream has more hits..
    auto collect = [&]() {
      if (trie_[i()].depth_and_hits.has_hit)
      {
        const auto needle_length = trie_[i()].depth_and_hits.depth;
        const auto text_start = text_pos - needle_length;
        // we want the first AAA of the spawn to be part of the hit; otherwise that hit will be reported by shorter sub-spawns or the master
        if (current_spawn_depths - needle_length >= spawn.max_prefix_loss_leftover) 
        {
          return false;
        }
        for (const auto needle_idx : umap_index2needles_.at(i()))
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
    if (i() == 0)
    {
      return 0;
    }
      
    // follow from the suffix...
    Index suf = trie_[i()].suffix;
    assert(suf.isValid());
    return follow_(suf, aa);
  }

  bool ACTrie::followSpawn_(ACSpawn& spawn, const AA edge, ACTrieState& state) const
  {
    // let spawn follow the original edge
    Index j = follow_(spawn.tree_pos, edge);
    const int new_depth = int(trie_[j()].depth_and_hits.depth);
    // did we loose a prefix?              old-depth                         new depth
    const int up_count = int(trie_[spawn.tree_pos()].depth_and_hits.depth) - new_depth + 1;
    if (up_count >= spawn.max_prefix_loss_leftover)
    { // spawn is dead because it lost its AAA/MM
      return false;
    }
    // update the prefix length
    spawn.max_prefix_loss_leftover -= up_count;
    spawn.tree_pos = j;
    addHitsSpawn_(j, spawn, spawn.textPos(state), state.hits, new_depth);
    return true;
  }

  Index ACTrie::stepMaster_(const Index i, const AA edge, ACTrieState& state) const
  {
    // has direct child (could also be an ambiguous AA - we don't care as long as a needle did contain that character)
    Index ch = findChildBFS_(i, edge);
    
    const bool consider_ambAA = max_aaa_ != 0;
    const bool consider_MM = max_mm_ != 0;
    
    // AAA
    if (edge.isAmbiguous())
    { // create spawns?
      AA from(edge), to(edge);
      if (consider_ambAA)
      { // first try AAA's (since edge is AAA)
        std::tie(from, to) = _getSpawnRange(edge);
        createSpawns_(i, from, to, state, max_aaa_ - 1, max_mm_);
      }
      // test all other AA's for mismatch
      if (consider_MM)
      { 
        createMMSpawns_(i, from, to, edge, state, max_aaa_, max_mm_ - 1);
      }

      // outdated dogma:
      // reset master to root (we are not allowed to have a path which contains an AAA in the master)
      //return {0};
    }                       
    // edge is unambiguous
    else if (consider_MM)
    { // try a MM for all AA's other than 'edge'
      createMMSpawns_(i, edge, edge, edge, state, max_aaa_, max_mm_ - 1);
    }
  
    // Master continues with the AA, no matter what it was...
    
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

  bool ACTrie::stepSpawn_(ACSpawn& spawn, ACTrieState& state) const
  {
    for (AA edge = spawn.nextValidAA(state); edge.isValid(); edge = spawn.nextValidAA(state))
    {
      const bool consider_ambAA = spawn.max_aaa_leftover > 0;
      const bool consider_MM = spawn.max_mm_leftover > 0;

      // AAA
      if (edge.isAmbiguous())
      { // create spawns from this spawn?
        AA from(edge), to(edge);
        if (consider_ambAA)
        { // first try AAA's (since edge is AAA)
          std::tie(from, to) = _getSpawnRange(edge);
          ACSpawn sp_temp = spawn;
          --sp_temp.max_aaa_leftover;
          createSubSpawns_(sp_temp, from, to, state);
        }
        // test all other superAA's for mismatch (except for AAA range, and the original edge itself)
        if (consider_MM)
        {
          ACSpawn sp_temp = spawn;
          --sp_temp.max_mm_leftover;
          createMMSubSpawns_(sp_temp, from, to, edge, state);
        }
      }
      else if (consider_MM) // edge is unambiguous
      { // try a MM for all superAA's other than 'edge'
        ACSpawn sp_temp = spawn;
        --sp_temp.max_mm_leftover;
        createMMSubSpawns_(sp_temp, edge, edge, edge, state);
      }

      // process the spawn itself
      if (!followSpawn_(spawn, edge, state)) return false;
    }
    return false; // end of query reached
  }

  void ACTrie::createMMSpawns_(const Index i, const AA except_fromAA, const AA except_toAA, const AA except_edge, ACTrieState& state, const uint32_t aaa_left, const uint32_t mm_left) const
  {
    // create super-AA range, i.e. including the ambiguous AA's, since a peptide could contain an 'X', which we would like to match
    auto [from, to] = _getSpawnRange(AA('$'));
    for (AA mm_aa = from; mm_aa <= to; ++mm_aa)
    {
      if (mm_aa == except_fromAA)
      { // ignore this range
        mm_aa = except_toAA;
        continue;
      }
      // ignore edge from spawn
      if (mm_aa == except_edge) 
      { 
        continue;
      }
      createSpawns_(i, mm_aa, mm_aa, state, aaa_left, mm_left);
    }
  }

  void ACTrie::createMMSubSpawns_(const ACSpawn& prototype, const AA except_fromAA, const AA except_toAA, const AA except_edge, ACTrieState& state) const
  {
    // create super-AA range, i.e. including the ambiguous AA's, since a peptide could contain an 'X', which we would like to match
    auto [from, to] = _getSpawnRange(AA('$'));
    for (AA mm_aa = from; mm_aa <= to; ++mm_aa)
    {
      if (mm_aa == except_fromAA)
      { // ignore this range
        mm_aa = except_toAA;
        continue;
      }
      // ignore edge from spawn
      if (mm_aa == except_edge)
      {
        continue;
      }
      createSubSpawns_(prototype, mm_aa, mm_aa, state);
    }
  }

  void ACTrie::createSpawns_(const Index i, const AA fromAA, const AA toAA, ACTrieState& state, const uint32_t aaa_left, const uint32_t mm_left) const
  {
    for (AA aa = fromAA; aa <= toAA; ++aa)
    {
      Index spawn_pos = follow_(i, aa); // call this using naive follow_(), which matches the exact char
      if (spawn_pos() > 0) // not at root
      {
        const uint8_t depth = trie_[spawn_pos()].depth_and_hits.depth;
        auto new_spawn = state.spawns.emplace(state.textPosIt(), // the master already points to the next AA, so spawn can start there
                                               spawn_pos,
                                               aaa_left,
                                               mm_left,
                                               depth);
        // we might have found a hit already: report it
        addHits_(spawn_pos, new_spawn.textPos(state), state.hits);
      }
    }
  }

  void ACTrie::createSubSpawns_(const ACSpawn& prototype, const AA fromAA, const AA toAA, ACTrieState& state) const
  {
    for (AA aa = fromAA; aa <= toAA; ++aa)
    {
      ACSpawn s(prototype);
      if (followSpawn_(s, aa, state))
      { // spawn survived following the edge
        state.spawns.push(std::move(s));
      }
    }
  }

  Index ACTrie::findChildNaive_(Index parent, AA child_label)
  {
    for (auto child : umap_index2children_naive_[parent])  // OMS_CODING_TEST_EXCLUDE Note: only a 4byte type. Copy it!
    {
      if (trie_[child.pos()].edge == child_label)
        return child;
    }
    return Index {};
  }

  Index ACTrie::findChildBFS_(const Index parent, const AA child_label) const
  {
    size_t start = trie_[parent()].first_child();
    size_t end = start + trie_[parent()].nr_children;
    for (size_t i = start; i < end; ++i)
    {
      if (trie_[i].edge == child_label)
        return Index::T(i);
    }
    return Index {};
  }

  void ACTrieState::setQuery(const std::string& haystack)
  {
    hits.clear();
    query_ = haystack;
    it_q_ = query_.begin();
    tree_pos = 0;
    while (!spawns.empty())
    {
      spawns.pop();
    }
  }

  size_t ACTrieState::textPos() const
  {
    return std::distance(query_.cbegin(), it_q_);
  }

  std::string::const_iterator ACTrieState::textPosIt() const
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
    return OpenMS::nextValidAA(query_.cend(), it_q_);
  }

  ACSpawn::ACSpawn(std::string::const_iterator query_pos, Index tree_pos, uint8_t max_aa, uint8_t max_mm, uint8_t max_prefix_loss) :
      it_query(query_pos), tree_pos(tree_pos), max_aaa_leftover(max_aa), max_mm_leftover(max_mm), max_prefix_loss_leftover(max_prefix_loss)
  {
  }

  size_t ACSpawn::textPos(const ACTrieState& state) const
  {
    return std::distance(state.query_.cbegin(), it_query);
  }

  AA ACSpawn::nextValidAA(const ACTrieState& state)
  {
    return OpenMS::nextValidAA(state.query_.cend(), it_query);
  }

  AA nextValidAA(const std::string::const_iterator end, std::string::const_iterator& it_q)
  {
    AA res {'?'}; // invalid
    while (it_q != end)
    {
      res = AA(*it_q);
      ++it_q;
      if (res.isValid())
      {
        return res;
      }
    }
    return res;
  }


  bool Index::isInvalid() const
  {
    return i_ == std::numeric_limits<T>::max();
  }

  bool Index::isValid() const
  {
    return i_ != std::numeric_limits<T>::max();
  }

  Index::T Index::operator()() const
  {
    return i_;
  }

  bool Index::operator==(const Index other) const
  {
    return i_ == other.i_;
  }

  Index::T& Index::pos()
  {
    return i_;
  }

  Index::T Index::pos() const
  {
    return i_;
  }

} // namespace OpenMS
