// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>
#include <bitset>
#include <cassert>
#include <cstring>
#include <functional> // for std::hash
#include <limits>
#include <queue>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <bit> // for std::popcount

namespace OpenMS
{
  /// representation of an AminoAcid (see AA class).
  /// Ambiguous AA's are consecutive (which saves effort during their enumeration)
  constexpr char const AAtoChar[28] = {
    'A', // 00 Ala Alanine
    'Y', // 01 Tyr Tyrosine
    'C', // 02 Cys Cysteine
    'D', // 03 Asp Aspartic Acid   // B
    'N', // 04 Asn Asparagine      // B
    'F', // 05 Phe Phenylalanine
    'G', // 06 Gly Glycine
    'H', // 07 His Histidine
    'I', // 08 Ile Isoleucine      // J
    'L', // 09 Leu Leucine         // J
    'K', // 10 Lys Lysine
    'W', // 11 Trp Tryptophan
    'M', // 12 Met Methionine
    'O', // 13 Pyl Pyrrolysine
    'P', // 14 Pro Proline
    'E', // 15 Glu Glutamic Acid   // Z
    'Q', // 16 Gln Glutamine       // Z
    'R', // 17 Arg Arginine
    'S', // 18 Ser Serine
    'T', // 19 Thr Threonine
    'U', // 20 Selenocysteine
    'V', // 21 Val Valine
    // ambiguous AAs start here (index: 22...25)
    'B', // 22 Aspartic Acid, Asparagine  $   // the ambAA's need to be consecutive (B,J,Z,X,$)
    'J', // 23 Leucine, Isoleucine        $
    'Z', // 24 Glutamic Acid, Glutamine   $
    'X', // 25 Unknown 
    // non-regular AA's, which are special
    '$', // 26 superAA, i.e. it models a mismatch, which can be anything, including AAAs
    '?', // 27 invalid AA (will usually be skipped) -- must be the last AA (AA::operator++ and others rely on it)
  };

  /// Conversion table from 7-bit ASCII char to internal value representation for an amino acid (AA)
  constexpr char const CharToAA[128] = {
    // ASCII char (7-bit Int with values from 0..127) --> amino acid 
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 0
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 1
  //                 $
    27, 27, 27, 27, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 2
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 3

  //  ,  A,  B,  C,  D,  E,  F,  G,  H,  I,  J,  K,  L,  M,  N,  O,
    27, 00, 22, 02, 03, 15, 05, 06, 07,  8, 23, 10,  9, 12, 04, 13, // 4

  // P,  Q,  R,  S,  T,  U,  V,  W,  X,  Y,  Z,   ,   ,   ,   ,   ,
    14, 16, 17, 18, 19, 20, 21, 11, 25, 01, 24, 27, 27, 27, 27, 27, // 5

  //  ,  a,  b,  c,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,
    27, 00, 22, 02, 03, 15, 05, 06, 07,  8, 23, 10,  9, 12, 04, 13,   // 6

  // p,  q,  r,  s,  t,  u,  v,  w,  x,  y,  z,   ,   ,   ,   ,   ,
    14, 16, 17, 18, 19, 20, 21, 11, 25, 01, 24, 27, 27, 27, 27, 27, // 7
  };

  /// Represents a needle found in the query.
  /// A needle (at position @p needle_index, as passed into ACTrie's addNeedle()) of length @p needle_length was found in the haystack (query) at position @p query_pos
  struct OPENMS_DLLAPI Hit {
    using T = uint32_t;
    Hit() = default;
    Hit(T needle_index, T needle_length, T query_pos) : needle_index(needle_index), needle_length(needle_length), query_pos(query_pos) {};
    T needle_index;
    T needle_length;
    T query_pos;
    bool operator==(const Hit& rhs) const
    {
      return needle_index == rhs.needle_index && needle_length == rhs.needle_length && query_pos == rhs.query_pos;
    }
  };

  /// An AminoAcid, which supports construction from char and has convenience functions such as
  /// isAmbiguous() or isValid()
  struct OPENMS_DLLAPI AA
  {
    /// Default C'tor; creates an invalid AA
    constexpr AA() : aa_(AA('?').aa_)
    {
    }
    
    /// C'tor from char; any char A-Z or a-z yields a valid AA.
    /// '$' is a special AA, which should only be used when modeling mismatches.
    /// All other chars produce an invalid AA ('?')
    constexpr explicit AA(const char c) : aa_(CharToAA[(unsigned char)c])
    {
    }

    /// yields the internal 8bit representation
    constexpr uint8_t operator()() const
    {
      return aa_;
    }

    /// equality operator
    constexpr bool operator==(const AA other) const
    {
      return aa_ == other.aa_;
    }

    /// less-or-equal operator
    constexpr bool operator<=(const AA other) const
    {
      return aa_ <= other.aa_;
    }

    /// is AA a 'B', 'J', 'Z', 'X', or '$' ?
    constexpr bool isAmbiguous() const
    {
      return aa_ >= AA('B').aa_;
    }
    
    /// is AA a letter or '$' ?
    constexpr bool isValid() const
    {
      return aa_ != AA('?').aa_;
    }

    /// is the AA a letter, i.e. A-Z or a-z?
    constexpr bool isValidForPeptide() const
    {
      return aa_ <= AA('X').aa_; // '$' or '?'
    }

    /// Pre-increment operator (advance to next AA)
    constexpr AA& operator++()
    {
      ++aa_;
      assert(aa_ <= AA('?').aa_); // make sure we don't overflow
      return *this;
    }

    /// Post-increment operator (advance to next AA)
    constexpr AA operator++(int)
    {
      AA r(*this);
      ++aa_;
      assert(aa_ <= AA('?').aa_); // make sure we don't overflow
      return r;
    }

    /// Decrement operator
    constexpr AA operator-(const AA rhs) const
    {
      AA r(*this);
      r.aa_ -= rhs.aa_;
      return r;
    }

    /// Convert this AA to a single letter representation, i.e. 'A'..'V', 'B', 'J', 'Z', 'X', '$', or '?' (invalid AA)
    constexpr char toChar() const
    {
      return AAtoChar[aa_];
    }

    /// Create an unambiguous AA from an index (0..21 inclusive -> 'A'..'V').
    /// @throw OpenMS::Precondition if not index < unambiguousAACount().
    static AA fromIndex(size_t index) noexcept
    {
      OPENMS_PRECONDITION(index < unambiguousAACount(), "AA::fromIndex(): index must be in [0, unambiguousAACount())");
      AA r;
      r.aa_ = static_cast<uint8_t>(index);
      return r;
    }

    /// Returns the number of unambiguous amino acids (A-Z, excluding B, J, Z, X), i.e. 22
    static size_t unambiguousAACount()
    {
      static_assert(AA('V').aa_ - AA('A').aa_ + 1 == 22, "Unambiguous AA count must be 22 (A-Z, excluding B, J, Z, X)");
      static_assert(AA('A').aa_ == 0, "A must be the first AA in the enumeration");
      static_assert(AA('V').aa_ == 21, "V must be the last unambiguous AA in the enumeration");
      return AA('V').aa_ - AA('A').aa_ + 1;
    }

  private:
    uint8_t aa_; ///< internal representation as 1-byte integer
  };

  /// An index with 32-bit representing the location of a node
  /// Allows to model invalid indices, see isInvalid() and isValid().
  class OPENMS_DLLAPI Index
  {
  public:
    using T = uint32_t;
    /// default C'tor; creates an invalid index
    constexpr Index() = default;
    
    /// C'tor from T
    constexpr Index(T val) : i_(val) {};

    /// is this Index invalid, i.e. should not be dereferenced
    constexpr bool isInvalid() const
    {
      return i_ == std::numeric_limits<T>::max();
    }

    /// is this Index valid, i.e. an actual index into a vector?
    constexpr bool isValid() const
    {
      return i_ != std::numeric_limits<T>::max();
    }

    /// convert to a number (might be invalid, check with .isValid() first)
    constexpr T operator()() const
    {
      return i_;
    }

    /// equality operator
    constexpr bool operator==(const Index other) const
    {
      return i_ == other.i_;
    }

    /// allows to set the index, using `index.pos() = 3;` or simply read its value
    constexpr T& pos()
    {
      return i_;
    }

    /// allows to read the index, using `index.pos()` 
    constexpr T pos() const
    {
      return i_;
    }
  private:
    T i_ = std::numeric_limits<T>::max(); ///< internal number representation; invalid state by default
  };
 } // namespace OpenMS

// this needs to go into namespace std
template<> struct std::hash<OpenMS::Index>
{
  std::size_t operator()(OpenMS::Index const& s) const noexcept
  {
    return std::hash<OpenMS::Index::T> {}(s());
  }
};

namespace OpenMS
{
  /// Custom bitset which uses a 32-bit integer to store bits (instead of 8 bytes for std::bitset<32> on Clang/GCC).
  struct Bitset {
    uint32_t bits = 0;

    // Set bit at position `pos` (0-based)
    void set(int pos) {
      bits |= (1u << pos);
    }

    // Clear bit at position `pos`
    void clear(int pos) {
      bits &= ~(1u << pos);
    }

    // Check if bit at position `pos` is set
    bool test(int pos) const {
      return (bits >> pos) & 1u;
    }

    // Left shift the bitset by `n` positions
    Bitset operator<<(int n) const {
      Bitset o = *this;
      return o <<= n;
    }

    // Left shift the bitset by `n` positions
    Bitset& operator<<=(int n) {
      if (n > 31)
      { // if we shift more than 31 bits, behaviour is undefined (e.g. GCC/Clang/MSVC will not modify 'bits' if n==32)
        bits = 0; // shift out all bits
      }
      else
      {
        bits <<= n;
      }
      return *this;
    }

    // Count number of set bits (population count)
    int pop_count() const {
      return std::popcount(bits);
    }
  }; // end Bitset


  /// A node in the AhoCorasick trie.
  /// Internally manages the suffix link and an index where its children start (this relies on the trie being stored in BFS order)
  struct OPENMS_DLLAPI ACNode
  {
    /// Default C'tor
    ACNode() {};

    /// C'tor from an edge @p label (from parent to this node) and a @p depth in the tree
    ACNode(const AA label, const uint8_t depth) : edge(label)
    {
      depth_and_hits.depth = depth;
    }

    /// internal struct to steal one bit from @p depth to use as hit indicator
    struct DepthHits {
      DepthHits()
      {
        memset(this, 0, sizeof *this); // make sure bitfields are 0; C++20 allows {} initialization ...
      };
      uint8_t has_hit: 1; ///< does a pattern end here (or when following suffix links)?
      // we could add another bit here to distinguish between a local hit and suffix hit, but on Windows, this slows it down
      uint8_t depth : 7; ///< depth of node in the trie
    };

    using ChildCountType = uint8_t; 

    // do not use std::bitset, since GCC/clang wastes 8 bytes here.
    //std::bitset<26> children_bitset; ///< bitfield of children (if tree is in BFS order); 26 bits are enough to cover all AA incl. (B,J,X,Z)
    Bitset children_bitset;  ///< bitfield of children (if tree is in BFS order); 26 bits are enough to cover all AA incl. (B,J,X,Z)
    Index suffix {0};        ///< which node is our suffix?
    Index first_child {0};   ///< which node contains our first child node (if tree is in BFS order)
    AA edge {0};             ///< what is the edge label (from parent to this node)
    ChildCountType nr_children = 0; ///< number of children (if tree is in BFS order); (popcount on children_bitset will do the same but we have the space due to padding)
    DepthHits depth_and_hits; ///< depth of node in the tree and one bit if a needle ends in this node or any of its suffices
  };
  
  // forward declaration
  struct ACTrieState;

  /// a spin-off search path through the trie, which can deal with ambiguous AAs and mismatches
  struct OPENMS_DLLAPI ACScout
  {
    /// No default C'tor
    ACScout() = delete;

    /// C'tor with arguments
    ACScout(const char* query_pos, Index tree_pos, uint8_t max_aa, uint8_t max_mm, uint8_t max_prefix_loss);

    /// Where in the text are we currently?
    size_t textPos(const ACTrieState& state) const;

    /// Return the next valid AA in the query. If the query was fully traversed, an invalid AA is returned.
    /// This moves the internal iterator for the query forwards.
    AA nextValidAA();

    const char* it_query = 0; ///< position in query
    Index tree_pos;                       ///< position in trie
    uint8_t max_aaa_leftover {0};         ///< number of ambiguous AAs the scout can yet tolerate before exceeding the limit
    uint8_t max_mm_leftover {0};          ///< number of mismatches the scout can yet tolerate before exceeding the limit
    uint8_t max_prefix_loss_leftover {0}; ///< number of AA's which can get lost by following suffix links, before the scout must retire; reaching 0 means retire
  };

  /// Return the first valid AA from current position of @p it_q in the string, or (if the string ends) an invalid AA
  /// On return, @p it_q points to AA after the returned AA.
  AA nextValidAA(const char*& it_q);

  /// A state object for an ACTrie, i.e. dynamic information when traversing the trie (which is 'const' after construction)
  /// Useful when using multi-threading; each thread can walk the trie and keep track of its state using an instance of this class
  struct OPENMS_DLLAPI ACTrieState
  {
    friend ACScout;
    /// Set a haystack (query) where the needles (patterns) are to be searched
    /// This also resets the current trie-node to ROOT, and voids the hits
    void setQuery(const std::string& haystack);

    /// Where in the text are we currently?
    size_t textPos() const;

    /// Where in the text are we currently?
    const char* textPosIt() const;

    /// The current query
    const std::string& getQuery() const;

    /// Return the next valid AA in the query. If the query was fully traversed, an invalid AA is returned.
    /// This moves the internal iterator for the query forwards.
    AA nextValidAA();

    std::vector<Hit> hits;             ///< current hits found
    std::queue<ACScout> scouts;        ///< initial scout points which are currently active and need processing
    Index tree_pos;                    ///< position in trie (for the Primary)
  private:
    const char* it_q_;                 ///< position in query
    std::string query_;                ///< current query ( = haystack = text)
  };

  /// An Aho Corasick trie (a set of nodes with suffix links mainly)
  class OPENMS_DLLAPI ACTrie
  {
  public:
    /** 
      @brief Default C'tor which just creates a root node

      @param max_aaa Maximum number of ambiguous amino acids (B,J,Z,X) allowed in a hit
      @param max_mm Maximum number of mismatched amino acids allowed in a hit
    */
    ACTrie(uint32_t max_aaa = 0, uint32_t max_mm = 0);

    /// D'tor
    ~ACTrie();

    /// Add a needle to build up the trie.
    /// Call compressTrie() after the last needle was added before searching
    /// @throw Exception::InvalidValue if @p needle contains an invalid amino acid (such as '*')
    void addNeedle(const std::string& needle);

    /// Convenience function; adds needles to build up the trie.
    /// Call compressTrie() after the last needle was added before searching
    /// @throw Exception::InvalidValue if @p needles contains an invalid amino acid (such as '*'); you can use getNeedleCount() to trace which needle did cause the exception
    void addNeedles(const std::vector<std::string>& needles);

    /// Convenience function which adds needles and immediately compresses the trie (i.e. no more needles can be added)
    /// @throw Exception::InvalidValue if @p needles contains an invalid amino acid (such as '*'); you can use getNeedleCount() to trace which needle did cause the exception
    void addNeedlesAndCompress(const std::vector<std::string>& needles);

    size_t size() const
    {
      return trie_.size();
    }

    /**
       @brief Traverses the trie in BFS order and makes it more compact and efficient to traverse

       Also creates the suffix links.

       Call this function after adding all needles, and before searching any queries.
    */
    void compressTrie();

    /// How many needles were added to the trie?
    size_t getNeedleCount() const;

    /// Set maximum number of ambiguous amino acids allowed during search.
    /// This must not be called in the middle of a search. Otherwise search results will be mixed.
    void setMaxAAACount(const uint32_t max_aaa);

    /// Maximum number of ambiguous amino acids allowed during search
    uint32_t getMaxAAACount() const;
    
    /// Set maximum number of mismatches allowed during search.
    /// This must not be called in the middle of a search. Otherwise search results will be mixed.
    void setMaxMMCount(const uint32_t max_mm);

    /// Maximum number of mismatches allowed during search
    uint32_t getMaxMMCount() const;

    /// Resume search at the last position in the query and node in the trie.
    /// If a node (or any suffices) are a hit, then @p state.hits is cleared & filled and true is returned.
    /// If the query ends and there is no hit, false is returned.
    bool nextHits(ACTrieState& state) const;

    /// Collects all hits from the current position in the query until the end of the query
    /// I.e. similar to while(next(state)) merge(hits_all, state.hits);
    void getAllHits(ACTrieState& state) const;

  private:
    /// Resume search at the last position in the query and node in the trie.
    /// If a node (or any suffices) are a hit, then @p state.hits is NOT cleared, but filled and true is returned.
    /// If the query ends and all scouts are processed, false is returned (but hits might still have changed)
    bool nextHitsNoClear_(ACTrieState& state) const;

    /// Insert a new child node into the trie (unless already present) when starting at parent node @p from and following the edge labeled @p edge.
    /// Return the index of the (new) child node.
    /// Note: This operates on the naive trie, not the BFS.
    Index add_(const Index from, const AA edge);

    /**
      @brief Add all hits occurring in node @p i (including all its suffix hits)

      @param i The ACNode where a needle ends (also all its suffices are checked)
      @param text_pos current position in query (i.e. end of matched hit)
      @param[out] hits Result vector which will be expanded with hits (if any)
      @return true if hits were found
    **/
    bool addHits_(Index i, const size_t text_pos, std::vector<Hit>& hits) const;

    /// same as addHits_, but only follows the suffix chain until the scout loses its prefix
    bool addHitsScout_(Index i, const ACScout& scout, const size_t text_pos, std::vector<Hit>& hits, const int current_scout_depths) const;

    /// Starting at node @p i, find the child with label @p edge
    /// This requires an exact match, e.g. if @p edge is an 'X' it will only match to 'X' in the trie (=needles)
    /// If no child exists, follow the suffix link and try again (until root is reached)
    /// Note: This operates on the BFS trie (after compressTrie()), not the naive one.
    Index follow_(const Index i, const AA edge) const;

    /// Advances @p scout by consuming @p edge; same as follow_(), just for a scout
    /// Returns true if scout is still alive, false otherwise
    bool followScout_(ACScout& scout, const AA edge, ACTrieState& state) const;

    /// Same as follow_, but considers trying mismatches and AAA's if possible (by adding scouts to @p state)
    /// @return The new tree node, where Primary is at after consuming @p edge
    Index stepPrimary_(const Index i, const AA edge, ACTrieState& state) const;

    /// Same as follow_, but considers trying mismatches and AAA's if possible (by adding scouts to @p state)
    /// @return true if scout is still alive, false otherwise (i.e. did loose its prefix)
    bool stepScout_(ACScout& scout, ACTrieState& state) const;

    /// Create scouts from an AAA or MM, starting at trie node @p i, following edges in range @p fromAA to @p toAA
    /// The number of AAA's/MM's left for the scout must be passed (these numbers already reflect the original edge label)
    void createScouts_(const Index i, const AA fromAA, const AA toAA, ACTrieState& state, const uint32_t aaa_left, const uint32_t mm_left) const;

    /// Create scouts from a scout with an AAA or MM, using @p prototype as template, following edges in range @p fromAA to @p toAA
    void createSubScouts_(const ACScout& prototype, const AA fromAA, const AA toAA, ACTrieState& state) const;

    /// Same as createScouts_, but instantiate all possible AA's except for those in the range from @p except_fromAA to @p except_toAA and the @p except_edge itself.
    void createMMScouts_(const Index i, const AA except_fromAA, const AA except_toAA, const AA except_edge, ACTrieState& state, const uint32_t aaa_left, const uint32_t mm_left) const;

    /// Same as createSubScouts_, but instantiate all possible AA's except for those in the range from @p except_fromAA to @p except_toAA and the @p except_edge itself.
    void createMMSubScouts_(const ACScout& prototype, const AA except_fromAA, const AA except_toAA, const AA except_edge, ACTrieState& state) const;

    /// During needle addition (naive trie), obtain the child with edge @p child_label from @p parent; if it does not exist, an invalid Index is returned
    Index findChildNaive_(Index parent, AA child_label);

    /// After compression (BFS trie), obtain the child with edge @p child_label from @p parent; if it does not exist, an invalid Index is returned
    Index findChildBFS_(const Index parent, const AA child_label) const;

    std::vector<ACNode> trie_;  ///< the trie, in either naive structure or BFS order (after compressTrie)
    uint32_t needle_count_ {0}; ///< total number of needles in the trie
    uint32_t max_aaa_ {0};      ///< maximum number of ambAAs allowed
    uint32_t max_mm_ {0};       ///< maximum number of mismatches allowed

    /// maps a node to which needles end there (valid for both naive and BFS tree)
    std::vector<std::vector<uint32_t>> vec_index2needles_ = { {} }; // one empty element to allow capacity doubling; note: sparse vector but still much faster and less memory than unordered_map
    /// maps the child nodes of a node for the naive tree; only needed during naive trie construction; storing children in the BFS trie modeled in the ACNodes directly
    std::vector<std::vector<Index>> vec_index2children_naive_ = { {} }; // one empty element to allow capacity doubling
  };

} // namespace OpenMS

