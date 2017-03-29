// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Sandro Andreotti, Tobias Rausch $
// --------------------------------------------------------------------------


#ifndef SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H
#define SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H

#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>

#ifdef NDEBUG
#define DEBUG_ONLY if (false)
#else
#define DEBUG_ONLY if (true)
#endif

namespace seqan
{

  //////////////////////////////////////////////////////////////////////////////
  // AhoCorasickAmbiguous
  //////////////////////////////////////////////////////////////////////////////

  /**
	@brief Extended Aho-Corasick algorithm capable of matching ambiguous amino acid in the pattern (i.e. proteins).

	...
	Features:
	  + blazingly fast
	  + low memory usage
	  + number of allowed ambAA's can be capped by user (default 3).


	This implementation is based on the original AC in SeqAn.

  */

  struct AhoCorasickAmbiguous_;
  typedef Tag<AhoCorasickAmbiguous_> AhoCorasickAmbiguous;

  //////////////////////////////////////////////////////////////////////////////

  /// state of an AC spawn, operating on a trie
  template <typename TNeedle>
  struct Spawn {
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TVertexDescriptor current_state;
    __uint8 max_DepthsDecrease; // maximum loss in depths of traversed nodes (both while reporting hits and changing its own state)

    Spawn() :
      current_state(getNil<TVertexDescriptor>()),
      max_DepthsDecrease(0)
    {}

    Spawn(TVertexDescriptor init_state, __uint8 current_depth) :
      current_state(init_state),
      max_DepthsDecrease(current_depth)
    {}
  };

  template <typename TNeedle>
  class Pattern<TNeedle, AhoCorasickAmbiguous> {
    //____________________________________________________________________________
  private:
    Pattern(Pattern const& other);
    Pattern const& operator=(Pattern const & other);
    //____________________________________________________________________________
  public:
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    // "constant" data, after construction of trie
	Holder<TNeedle> data_host; // holds needles
    String<String<TSize> > data_terminalStateMap; // regular trie data -- plus: this gets augmented with all suffix traversals which are output nodes
    TGraph data_graph;                            // regular trie data
    String<TVertexDescriptor> data_supplyMap;     // trie suffix links
    String<__uint8> data_nodeDepths;              // depths of each graph node

    // To restore the automaton after a hit
    String<TSize> data_endPositions;	// All remaining keyword indices
    TSize data_keywordIndex;			// Current keyword that produced a hit
    TSize data_needleLength;			// Last length of needle to reposition finder
    TVertexDescriptor data_lastState;   // Last state of master instance in the trie
    typedef typename std::list<Spawn<TNeedle> > Spawns;
    Spawns spawns; // AC instances currently walking the tree

    //____________________________________________________________________________
    Pattern() {}

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
      SEQAN_CHECKPOINT
      setHost(*this, ndl);
    }

    ~Pattern() {
      SEQAN_CHECKPOINT
    }
    //____________________________________________________________________________
  };


  //////////////////////////////////////////////////////////////////////////////
  // Host Metafunctions
  //////////////////////////////////////////////////////////////////////////////

  template <typename TNeedle>
  struct Host< Pattern<TNeedle, AhoCorasickAmbiguous> >
  {
    typedef TNeedle Type;
  };

  template <typename TNeedle>
  struct Host< Pattern<TNeedle, AhoCorasickAmbiguous> const>
  {
    typedef TNeedle const Type;
  };


  //////////////////////////////////////////////////////////////////////////////
  // Functions
  //////////////////////////////////////////////////////////////////////////////

  template <typename TNeedle>
  inline void _createAcTrie(Pattern<TNeedle, AhoCorasickAmbiguous> & me)
  {

    //OpenMS::StopWatch sw;
    //sw.start();

    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

    //foo(TPosition());
    //foo(TKeyword());
    //foo(TAlphabet());
    //foo(TGraph());

    // Create regular trie
    createTrie(me.data_graph, me.data_terminalStateMap, host(me));

    // Create parent map
    String<TVertexDescriptor> parentMap;  //< allows to find parent of each node
    String<TAlphabet> parentCharMap;      //< allows to find character that led us to current node
    resizeVertexMap(me.data_graph, parentMap);
    resizeVertexMap(me.data_graph, parentCharMap);
    // init all nodes to Nil
    for (TPosition i = 0; i < length(parentMap); ++i) {
      assignProperty(parentMap, i, nilVal);
    }
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itEd(me.data_graph);
    //int count_prop = 0;
    for (;!atEnd(itEd);goNext(itEd)) {
      //              property,  vertex ,  value
      assignProperty(parentMap, targetVertex(itEd), sourceVertex(itEd));
      assignProperty(parentCharMap, targetVertex(itEd), label(itEd));
      //++count_prop;
    }
    //std::cout << "\nassigned " << count_prop << " props to " << length(parentMap) << " values\n";


    // Build AC
    TVertexDescriptor root = getRoot(me.data_graph);
    // properties....
    resizeVertexMap(me.data_graph, me.data_supplyMap);  // suffix links
    assignProperty(me.data_supplyMap, root, nilVal);
    resizeVertexMap(me.data_graph, me.data_nodeDepths);  // node depths
    assignProperty(me.data_nodeDepths, root, 0);

    // Bfs Traversal
    typedef typename Iterator<TGraph, BfsIterator>::Type TBfsIterator;
    TBfsIterator it(me.data_graph, root);
    //std::map<int, int> connectivity;
    goNext(it); // skip root

    for (;!atEnd(it);goNext(it))
    {
      const typename GetValue<TBfsIterator>::Type itval = *it; // dereferencing *it will give an index into the property array!

      //int edgecount = outDegree(me.data_graph, itval);
      //++connectivity[edgecount];

      TVertexDescriptor parent = getProperty(parentMap, itval);
      assignProperty(me.data_nodeDepths, itval, getProperty(me.data_nodeDepths, parent) + 1);

      TAlphabet sigma = getProperty(parentCharMap, itval);
      TVertexDescriptor down = getProperty(me.data_supplyMap, parent);
      while ((down != nilVal) &&
        (getSuccessor(me.data_graph, down, sigma) == nilVal))
      {
        down = getProperty(me.data_supplyMap, down);
      }
      if (down != nilVal)
      {
        assignProperty(me.data_supplyMap, itval, getSuccessor(me.data_graph, down, sigma));
        String<TPosition> endPositions = getProperty(me.data_terminalStateMap, getProperty(me.data_supplyMap, itval));
        if (!empty(endPositions))
        {
          String<TPosition> endPositionsCurrent = getProperty(me.data_terminalStateMap, itval);
          typedef typename Iterator<String<TPosition>, Rooted >::Type TStringIterator;
          TStringIterator sit = begin(endPositions);
          for (;!atEnd(sit); goNext(sit))
          {
            appendValue(endPositionsCurrent, *sit);
          }
          assignProperty(me.data_terminalStateMap, itval, endPositionsCurrent);
        }
      }
      else {
        assignProperty(me.data_supplyMap, itval, root);
      }

    }

    //sw.stop();
    //std::cout << OpenMS::String(" Trie build time: ") + sw.getClockTime() + " s (wall), " + sw.getCPUTime() + " s (CPU)).\n\n";

    //for (std::map<int,int>::const_iterator it = connectivity.begin(); it != connectivity.end(); ++it)
    //{
    //  std::cout << "   conn: " << it->first << " : " << it->second << "x\n";
    //}

  }

  template <typename TNeedle, typename TNeedle2>
  void setHost(Pattern<TNeedle, AhoCorasickAmbiguous> & me, TNeedle2 const & needle) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT(empty(needle));
    setValue(me.data_host, needle);
    clear(me.data_graph);
    clear(me.data_supplyMap);
    clear(me.data_endPositions);
    clear(me.data_terminalStateMap);
    _createAcTrie(me);
    me.data_needleLength = 0;
    me.spawns.clear();

    //fstream strm;
    //strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
    //String<String<char> > nodeMap;
    //_createTrieNodeAttributes(me.data_graph, me.data_terminalStateMap, nodeMap);
    //String<String<char> > edgeMap;
    //_createEdgeAttributes(me.data_graph,edgeMap);
    //write(strm,me.data_graph,nodeMap,edgeMap,DotDrawing());
    //strm.close();
    // Supply links
    //for(unsigned int i=0;i<length(me.data_supplyMap);++i) {
    //	std::cout << i << "->" << getProperty(me.data_supplyMap,i) << ::std::endl;
    //}
  }
  template <typename TNeedle, typename TNeedle2>
  inline void setHost(Pattern<TNeedle, AhoCorasickAmbiguous> & me, TNeedle2 & needle)
  {
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
  }
  //____________________________________________________________________________
  template <typename TNeedle>
  inline void _patternInit(Pattern<TNeedle, AhoCorasickAmbiguous> & me)
  {
    SEQAN_CHECKPOINT
    clear(me.data_endPositions);
    me.data_keywordIndex = 0;
	me.data_needleLength = 0;
    me.data_lastState = getRoot(me.data_graph);
	me.spawns.clear();
  }
  //____________________________________________________________________________
  template <typename TNeedle>
  inline typename Size<TNeedle>::Type position(Pattern<TNeedle, AhoCorasickAmbiguous> & me)
  {
    return me.data_keywordIndex;
  }

  template <typename TFinder, typename TNeedle>
  inline void _reportHit(TFinder & finder, Pattern<TNeedle, AhoCorasickAmbiguous> & me) {
    size_t idx_endPosVec = length(me.data_endPositions) - 1;
    me.data_keywordIndex = me.data_endPositions[idx_endPosVec];
    resize(me.data_endPositions, idx_endPosVec); // pop last hit
    me.data_needleLength = length(value(host(me), me.data_keywordIndex)) - 1;
    finder -= me.data_needleLength; // position finder at beginning of hit
    _setFinderLength(finder, me.data_needleLength + 1);
    _setFinderEnd(finder, position(finder) + length(finder)); // end of match within haystack
    return;
  }

  template<typename T>
  inline bool isAmbiguous(T c)
  {
    static const typename ValueSize<T>::Type vB = ordValue(T('B')); // D,N
    static const typename ValueSize<T>::Type vZ = ordValue(T('Z')); // E,Q
    static const typename ValueSize<T>::Type vX = ordValue(T('X')); // all
    static const uint32_t anyAA = (1 << vB) | (1 << vZ) | (1 << vX);
    return ((1 << ordValue(c)) & anyAA);
  }

  /**
   @brief given an ambAA @p c, return a range of AA's which need to be spawned and an extra AA which is meant for the master thread.
  */
  template<typename T>
  inline AminoAcid _getSpawnRange(T& idxFirst, T& idxLast, const AminoAcid c)
  {
    // jump table:                 // AA for main thread     // start of spawns        // end of spawns (including)
    static const T jump[3][3] = { { ordValue(AminoAcid('N')), ordValue(AminoAcid('D')), ordValue(AminoAcid('D')) },  // B = D,N
                                  { ordValue(AminoAcid('Q')), ordValue(AminoAcid('E')), ordValue(AminoAcid('E')) },  // Z = E,Q
                                  { 0,                        1,                        ordValue(AminoAcid('V')) } };// X = A..V
    SEQAN_ASSERT(ordValue(AminoAcid('B')) == 20); // make sure the table is ordered as we expect
    idxFirst = jump[ordValue(c) - 20][1];
    idxLast = jump[ordValue(c) - 20][2];
    return AminoAcid(jump[ordValue(c) - 20][0]); // AA for main thread
  }

  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline AminoAcid _createSpawns(Pattern<TNeedle, AhoCorasickAmbiguous>& me,
                                 typename Pattern<TNeedle, AhoCorasickAmbiguous>::TVertexDescriptor& current,
                                 AminoAcid c)
  {
    if (isAmbiguous(c))
    { // spawn children
	  DEBUG_ONLY std::cout << "found AAA: " << c << "\n";
      typedef typename Size<AminoAcid>::Type TSize;
      TSize idxFirst, idxLast;
      c = _getSpawnRange(idxFirst, idxLast, c);
      for (TSize i = idxFirst; i <= idxLast; ++i)
      {
        typedef typename Pattern<TNeedle, AhoCorasickAmbiguous>::TVertexDescriptor TVertexDescriptor;
        TVertexDescriptor spawn_state = me.data_lastState; // a potential spawn
        if (_consumeChar(me, spawn_state, AminoAcid(i)))
        {
		  DEBUG_ONLY std::cout << "Spawn '" << AminoAcid(i) << "' created\n";
          __uint8 node_depth = getProperty(me.data_nodeDepths, spawn_state); // depth at which the AA was consumed!
          me.spawns.push_front(Spawn<TNeedle>(spawn_state, node_depth));
        }
      }
    }
    return c;
  }
  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline AminoAcid _createSpawns(Pattern<TNeedle, AhoCorasickAmbiguous>& me,
                                 const Spawn<TNeedle>& current,
                                 AminoAcid c)
  {
    if (isAmbiguous(c))
    { // spawn children
	  DEBUG_ONLY std::cout << "found AAA: " << c << "\n";
      typedef typename Size<AminoAcid>::Type TSize;
      TSize idxFirst, idxLast;
      c = _getSpawnRange(idxFirst, idxLast, c);
      for (TSize i = idxFirst; i <= idxLast; ++i)
      {
        Spawn<TNeedle> spawn = current; // a potential spawn
        if (_consumeChar(me, spawn, AminoAcid(i)))
        {
		  DEBUG_ONLY std::cout << "Spawn '" << AminoAcid(i) << "' created\n";
          //__uint8 node_depth = getProperty(me.data_nodeDepths, spawn_state); // depth at which the AA was consumed!
          // Spawns of Spawns inherit the depths from their parents, since the master will also see this same AAA and spawn himself
          me.spawns.push_front(spawn);
        }
      }
    }
    return c;
  }

  // returns false if it reached the 'root', true otherwise;
  // might create new spawns
  template <typename TNeedle>
  inline bool _consumeChar(Pattern<TNeedle, AhoCorasickAmbiguous>& me,
                           Spawn<TNeedle>& spawn,
                           AminoAcid c)
  {
    typedef typename Pattern<TNeedle, AhoCorasickAmbiguous>::TVertexDescriptor TVertexDescriptor;

    c = _createSpawns(me, spawn, c); // spawn children and leave fixed AA for master spawn

    TVertexDescriptor successor, suffix_node;
    static const TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    // if we cannot go down, and up is possible:
    while (((successor = getSuccessor(me.data_graph, spawn.current_state, c)) == nilVal) &&
           ((suffix_node = getProperty(me.data_supplyMap, spawn.current_state)) != nilVal))
    { // .. follow suffix links upwards (if we do not fall below depth limit)
      if (false && suffix_node == getRoot(me.data_graph)) // not strictly required, since max_depths would work as well
      {
		DEBUG_ONLY std::cout << "spawn died while trying to match '" << c << "' passing root while going up\n";
        return false; // this spawn just threw away its reason of existance (it has no history). Die!
      }
      __uint8 depthDiff = getProperty(me.data_nodeDepths, spawn.current_state) - getProperty(me.data_nodeDepths, suffix_node);
      if (spawn.max_DepthsDecrease <= depthDiff) {
		DEBUG_ONLY std::cout << "spawn died while trying to match '" << c << "' and going up (AAA out of scope)\n";
        return false; // this spawn just threw away its reason of existance (i.e. the AAA). Die!
      }
      spawn.max_DepthsDecrease -= depthDiff;
      spawn.current_state = suffix_node;
    }
    // found a successor:
    if (successor != nilVal) {
      spawn.current_state = successor;
	  DEBUG_ONLY std::cout << "spawn matched '" << c << "'\n";

      // get hits
      typedef typename Size<TNeedle>::Type TSize;
      String<TSize> needle_hits = getProperty(me.data_terminalStateMap, spawn.current_state);
      if (length(needle_hits) > 0)
      {
        int path_length = getProperty(me.data_nodeDepths, spawn.current_state); // == length of current path to spawn
        int unambiguous_suffix_length = path_length - spawn.max_DepthsDecrease; // == length of suffix peptide which does not contain AAA
		DEBUG_ONLY std::cout << "  spawn adding hits which are more than '" << unambiguous_suffix_length << "' chars long (thus contain the AAA).\n";

        // but only report those which contain the AAA
        for (int i = 0; i < length(needle_hits); ++i)
        {
          int hit_length = length(value(host(me), needle_hits[i]));
          if (hit_length <= unambiguous_suffix_length) break; // assumption: terminalStateMap is sorted by length of hits! ... uiuiui...
		  DEBUG_ONLY std::cout << "  spawn hit: #" << i << "\n"; // value(value(host(me), needle_hits[i]))
          append(me.data_endPositions, needle_hits[i]); // append hits which still contain the AAA
        }
      }
      
      return true;
    }
    else { // no fitting successor and no way upwards --> we are at root
      //spawn.current_state = getRoot(me.data_graph); // not required -- this Spawn will die anyway
      return false;
    }
  }

  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _consumeChar(Pattern<TNeedle, AhoCorasickAmbiguous>& me,
                           typename Pattern<TNeedle, AhoCorasickAmbiguous>::TVertexDescriptor& current,
                           AminoAcid c)
  {
    typedef typename Pattern<TNeedle, AhoCorasickAmbiguous>::TVertexDescriptor TVertexDescriptor;
    
    c = _createSpawns(me, current, c); // and leave fixed AA for master thread
    
    TVertexDescriptor successor;
    static const TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    // if we cannot go down, and up is possible:
    while (((successor = getSuccessor(me.data_graph, current, c)) == nilVal) &&
      (getProperty(me.data_supplyMap, current) != nilVal))
    { // .. follow suffix links upwards
      current = getProperty(me.data_supplyMap, current);
    }
    // found a successor:
    if (successor != nilVal) {
      current = successor;
	  DEBUG_ONLY std::cout << "main thread/spawntest consumed '" << c << "'\n";
      typedef typename Size<TNeedle>::Type TSize;
      String<TSize> needle_hits = getProperty(me.data_terminalStateMap, current);
      if (length(needle_hits))
      {
		DEBUG_ONLY std::cout << "  hit count: '" << length(needle_hits) << "'\n";
        append(me.data_endPositions, getProperty(me.data_terminalStateMap, current)); // indices into TNeedle!
      }
      return true;
    }
    else { // no fitting successor and no way upwards --> we are at root
	  DEBUG_ONLY std::cout << "main thread/spawntest remains at root after consuming '" << c << "'\n";
      current = getRoot(me.data_graph);
      return false;
    }
  }

  template <typename T = void>
  struct AAEquivalenceClass_
  {
    static unsigned const VALUE[24]; // EquivalenceClassAA_<>::VALUE
  };
  template <typename T>
  unsigned const AAEquivalenceClass_<T>::VALUE[24] =
  {
    1, // 0 Ala Alanine (A)
    2, // 1 Arg Arginine (R)
    4, // 2 Asn Asparagine (N)
    8, // 3 Asp Aspartic Acid (D)
    16, // 4 Cys Cystine (C)
    32, // 5 Gln Glutamine (Q)
    64, // 6 Glu Glutamic Acid (E)
    128, // 7 Gly Glycine (G)
    256, // 8 His Histidine (H)
    512, // 9 Ile Isoleucine (I)
    1024, // 10 Leu Leucine (L)
    2048, // 11 Lys Lysine (K)
    4096, // 12 Met Methionine (M)
    8192, // 13 Phe Phenylalanine (F)
    16384, // 14 Pro Proline (P)
    32768, // 15 Ser Serine (S)
    65536, // 16 Thr Threonine (T)
    131072, // 17 Trp Tryptophan (W)
    262144, // 18 Tyr Tyrosine (Y)
    524288, // 19 Val Valine (V)
    4 + 8, // 20  Aspartic Acid (D), Asparagine(N) == (B)
    32 + 64, // 21 Glutamic Acid(E), Glutamine(Q) == (Z)
    static_cast<unsigned>(2^(19+1) - 1), // 22 (X) matches ALL
    static_cast<unsigned>(-1), // 23 Terminator (dummy)
  };

  template <typename TFinder, typename TNeedle>
  inline bool find(TFinder & finder, Pattern<TNeedle, AhoCorasickAmbiguous> & me) {
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    if (empty(finder)) {
      _patternInit(me);
      _finderSetNonEmpty(finder);
    }
    else {
      finder += me.data_needleLength; // restore last consumed position in haystack
      if (!empty(me.data_endPositions)) { // Process left-over hits
        _reportHit(finder, me);
        return true;
      }
      ++finder; // advance to next position
    }

    while (!atEnd(finder)) {
      const AminoAcid c = *finder;
      // spawns; do them first, since we might add new spawns in main thread which are however settled at that point
      if (!me.spawns.empty()) {
        typedef typename Pattern<TNeedle, AhoCorasickAmbiguous>::Spawns::iterator SpawnCIt;
        SpawnCIt it = me.spawns.begin();
        while (it != me.spawns.end())
        {
          if (!_consumeChar(me, *it, c)) // might create new spawns
          { // spawn reached root --> kill it
            it = me.spawns.erase(it); // remove and advance
          }
          else ++it; // next spawn
        }
      }
      // main thread
      _consumeChar(me, me.data_lastState, c); // might create new spawns
      
      if (!empty(me.data_endPositions)) {
        _reportHit(finder, me);
        return true;
      }

      ++finder;
    }
    return false;
  }

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H
