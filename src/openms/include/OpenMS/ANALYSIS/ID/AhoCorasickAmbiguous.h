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
  // AhoCorasickAmb
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

  struct AhoCorasickAmb_;
  typedef Tag<AhoCorasickAmb_> AhoCorasickAmb;

  //////////////////////////////////////////////////////////////////////////////

  /// state of an AC spawn, operating on a trie
  template <typename TNeedle>
  struct Spawn {
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVert;
    typedef typename Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    TVert current_state;
    KeyWordLengthType max_DepthsDecrease; // maximum loss in depths of traversed nodes (both while reporting hits and changing its own state)
    KeyWordLengthType ambAA_seen;         // number of ambAA's which the spawn has seen

    private:
	    Spawn();

    public:
	    Spawn(TVert init_state, KeyWordLengthType current_depth, KeyWordLengthType aaa_seen) :
		    current_state(init_state),
		    max_DepthsDecrease(current_depth),
        ambAA_seen(aaa_seen)
	    {}	
  };

  template <typename TNeedle>
  class Pattern<TNeedle, AhoCorasickAmb> {
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
    typedef typename VertexDescriptor<TGraph>::Type TVert;
    typedef __uint8 KeyWordLengthType;


    // constant after C'tor; 
    const TVert nilVal;
    const KeyWordLengthType max_ambAA;  // default: 3
                                        
    // "constant" data, after construction of trie
	  Holder<TNeedle> data_host; // holds needles
    String<String<TSize> > data_terminalStateMap; // regular trie data -- plus: this gets augmented with all suffix traversals which are output nodes
    TGraph data_graph;                            // regular trie data
    String<TVert> data_supplyMap;     // trie suffix links
    String<KeyWordLengthType> data_nodeDepths;              // depths of each graph node

    // "working" set; changes with every hit
    String<TSize> data_endPositions;	// All remaining keyword indices
    TSize data_keywordIndex;			// Current keyword that produced a hit
    TSize data_needleLength;			// Last length of needle to reposition finder
    TVert data_lastState;   // Last state of master instance in the trie
    typedef typename std::list<Spawn<TNeedle> > Spawns;
    typedef typename std::list<Spawn<TNeedle> >::iterator SpawnIt;
    typedef typename std::list<Spawn<TNeedle> >::const_iterator SpawnCIt;
    Spawns spawns;                      // spawn instances currently walking the tree
    typedef typename std::list<KeyWordLengthType> AmbAAPositions;
    typedef typename std::list<KeyWordLengthType>::iterator AmbAAPositionsIt;
    AmbAAPositions ambAA_positions;    // indices of ambAA's relative to current path in trie; when going up, this list must be updated

    //____________________________________________________________________________
    Pattern() {}

    template <typename TNeedle>
    Pattern(TNeedle const & ndl, KeyWordLengthType max_AAA = 3)
      : nilVal(getNil<TVert>()),
        max_ambAA(max_AAA)
    {
      SEQAN_CHECKPOINT
      typedef typename Value<TNeedle>::Type TKeyword;
      typedef typename Value<TKeyword>::Type TAlphabet;
      
      LOG_INFO << "AC with " << int(max_ambAA) << " ambiguous AAs!" << std::endl;

      for (TSize i = 0; i < length(ndl); ++i)
      {
        if (length(ndl[i]) > numeric_limits<KeyWordLengthType>::max())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, "Pattern<AhoCorasickAmb>(PeptideSet)", std::string("Input peptide to AhoCorasickAmb must NOT be longer than 255 chars!").c_str(), std::string(begin(ndl[i]), end(ndl[i])));
        }
        for (TSize j = 0; j < length(ndl[i]); ++j)
        {
          if (isAmbiguous(ndl[i][j]))
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, "Pattern<AhoCorasickAmb>(PeptideSet)", std::string("Input peptide to AhoCorasickAmb must NOT contain ambiguous amino acids ('B'/'Z'/'X')! Note: unknown AAs (e.g. 'U') will be converted to 'X' implicitly!").c_str(), std::string(begin(ndl[i]), end(ndl[i])));
          }
        }
      }
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
  struct Host< Pattern<TNeedle, AhoCorasickAmb> >
  {
    typedef TNeedle Type;
  };

  template <typename TNeedle>
  struct Host< Pattern<TNeedle, AhoCorasickAmb> const>
  {
    typedef TNeedle const Type;
  };


  //////////////////////////////////////////////////////////////////////////////
  // Functions
  //////////////////////////////////////////////////////////////////////////////

  template <typename TNeedle>
  inline void _createAcTrie(Pattern<TNeedle, AhoCorasickAmb> & me)
  {
    //OpenMS::StopWatch sw;
    //sw.start();
    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVert;
    TVert nilVal = getNil<TVert>();

    // Create regular trie
    createTrie(me.data_graph, me.data_terminalStateMap, host(me));

    // Create parent map
    String<TVert> parentMap;  //< allows to find parent of each node
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
      //             property,  vertex            , value
      assignProperty(parentMap, targetVertex(itEd), sourceVertex(itEd));
      assignProperty(parentCharMap, targetVertex(itEd), label(itEd));
      //++count_prop;
    }
    //std::cout << "\nassigned " << count_prop << " props to " << length(parentMap) << " values\n";


    // Build AC
    TVert root = getRoot(me.data_graph);
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

      TVert parent = getProperty(parentMap, itval);
      assignProperty(me.data_nodeDepths, itval, getProperty(me.data_nodeDepths, parent) + 1);

      TAlphabet sigma = getProperty(parentCharMap, itval);
      TVert down = getProperty(me.data_supplyMap, parent);
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
  void setHost(Pattern<TNeedle, AhoCorasickAmb> & me, TNeedle2 const & needle) {
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
  inline void setHost(Pattern<TNeedle, AhoCorasickAmb> & me, TNeedle2 & needle)
  {
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
  }
  //____________________________________________________________________________
  template <typename TNeedle>
  inline void _patternInit(Pattern<TNeedle, AhoCorasickAmb> & me)
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
  inline typename Size<TNeedle>::Type position(Pattern<TNeedle, AhoCorasickAmb> & me)
  {
    return me.data_keywordIndex;
  }

  template <typename TFinder, typename TNeedle>
  inline void _reportHit(TFinder & finder, Pattern<TNeedle, AhoCorasickAmb> & me) {
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
    OPENMS_PRECONDITION(isAmbiguous(c), "AminoAcid is not ambiguous! InternalError!");

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
  // This will only be called by the master itself, and at this point, we can surely add another ambAA (goUp() was called before)!
  template <typename TNeedle>
  inline AminoAcid _createSpawns(Pattern<TNeedle, AhoCorasickAmb>& me,
                                 typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current,
                                 AminoAcid c) // ALWAYS ambiguous!!!!
  {
	  DEBUG_ONLY std::cout << "found AAA: " << c << "\n";
    typedef typename Size<AminoAcid>::Type TSize;
    typedef Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    TSize idxFirst, idxLast;
    c = _getSpawnRange(idxFirst, idxLast, c);
    for (TSize i = idxFirst; i <= idxLast; ++i)
    {
      typedef typename Pattern<TNeedle, AhoCorasickAmb>::TVert TVert;
      TVert node_spawn = me.data_lastState; // a potential spawn
      if (_consumeChar(me, node_spawn, AminoAcid(i), Tag<FixedAASpec>()))
      {
		    DEBUG_ONLY std::cout << "Spawn '" << AminoAcid(i) << "' created\n";
        const KeyWordLengthType node_depth = getProperty(me.data_nodeDepths, node_spawn); // depth at which the AA was consumed!
        // push_front is paramount, since we might iterate over old spawns at this very moment
        me.spawns.push_front(Spawn<TNeedle>(node_spawn, node_depth, (KeyWordLengthType)me.ambAA_positions.size()));
      }
    }
    return c;
  }
  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline AminoAcid _createSpawns(Pattern<TNeedle, AhoCorasickAmb>& me,
                                 const Spawn<TNeedle>& current,
                                 AminoAcid c) // ALWAYS ambiguous!!!!
  {
	  DEBUG_ONLY std::cout << "found AAA: " << c << "\n";
    typedef typename Size<AminoAcid>::Type TSize;
    typedef Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    TSize idxFirst, idxLast;
    c = _getSpawnRange(idxFirst, idxLast, c);
    for (TSize i = idxFirst; i <= idxLast; ++i)
    {
      Spawn<TNeedle> spawn = current; // a potential spawn
      if (_consumeChar(me, spawn, AminoAcid(i), Tag<FixedAASpec>()))
      {
		    DEBUG_ONLY std::cout << "Spawn '" << AminoAcid(i) << "' created\n";
        //KeyWordLengthType node_depth = getProperty(me.data_nodeDepths, spawn_state); // depth at which the AA was consumed!
        // Spawns of Spawns inherit the depths from their parents, since the master will also see this same AAA and spawn himself
        me.spawns.push_front(spawn);
      }
    }
    return c;
  }



  struct IsAmbAASpec {};
  struct FixedAASpec {};

  /// ###  go down  ####
  template<class TNeedle> inline bool goDown(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current, AminoAcid c)
  { 
    Pattern<TNeedle, AhoCorasickAmb>::TVert successor = getSuccessor(me.data_graph, current, c);
    if (successor == me.nilVal) return false;
    DEBUG_ONLY std::cout << "master/test matched '" << c << "'\n";
    current = successor;
    return true;
  }
  template<class TNeedle> inline bool goDown(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Spawn<TNeedle>& spawn, AminoAcid c)
  {
    Pattern<TNeedle, AhoCorasickAmb>::TVert successor = getSuccessor(me.data_graph, spawn.current_state, c);
    if (successor == me.nilVal) return false;
    DEBUG_ONLY std::cout << "spawn matched '" << c << "'\n";
    spawn.current_state = successor;
    return true;
  }

  /// ###  go up  ####
  template<class TNeedle> inline bool goUp(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Spawn<TNeedle>& spawn)
  {
    if (atRoot(me, spawn)) return false;
    Pattern<TNeedle, AhoCorasickAmb>::TVert suffix_node = getProperty(me.data_supplyMap, spawn.current_state);
    // check if spawn is allowed to loose that many chars in front
    typedef Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    KeyWordLengthType depthDiff = getProperty(me.data_nodeDepths, spawn.current_state) - getProperty(me.data_nodeDepths, suffix_node);
    if (spawn.max_DepthsDecrease <= depthDiff) {
      DEBUG_ONLY std::cout << "spawn died while going up (AAA out of scope)\n";
      spawn.current_state = getRoot(me.data_graph); // reset to root -- indicating failure!
      return false; // this spawn just threw away its reason of existance (i.e. the AAA). Die!
    }
    spawn.max_DepthsDecrease -= depthDiff;
    spawn.current_state = suffix_node; // no need to check for nilVal, since we cannot reach root (depths runs out before!)
    return true;
  }
  template<class TNeedle> inline bool goUp(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state)
  {
    if (atRoot(me, current_state)) return false;
    Pattern<TNeedle, AhoCorasickAmb>::TVert suffix_node = getProperty(me.data_supplyMap, current_state);
    if (suffix_node != me.nilVal) {
      current_state = suffix_node;
      return true;
    }
    return false;
  }

  template<class TNeedle> inline bool atRoot(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Spawn<TNeedle>& spawn) {return spawn.current_state == getRoot(me.data_graph);}
  template<class TNeedle> inline bool atRoot(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state)
  {
    DEBUG_ONLY std::cout << "master/test remains at root\n";
    return current_state == getRoot(me.data_graph);
  }

  template<class TNeedle> inline void addHits(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Spawn<TNeedle>& spawn)
  {
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
        int hit_length = (int)length(value(host(me), needle_hits[i]));
        if (hit_length <= unambiguous_suffix_length) break; // assumption: terminalStateMap is sorted by length of hits! ... uiuiui...
        DEBUG_ONLY std::cout << "  spawn hit: #" << i << "\n"; // value(value(host(me), needle_hits[i]))
        append(me.data_endPositions, needle_hits[i]); // append hits which still contain the AAA
      }
    }
  }
  template<class TNeedle> inline void addHits(typename Pattern<TNeedle, AhoCorasickAmb>& me, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state)
  {
    typedef typename Size<TNeedle>::Type TSize;
    String<TSize> needle_hits = getProperty(me.data_terminalStateMap, current_state);
    if (length(needle_hits))
    {
      DEBUG_ONLY std::cout << "master/test hit total count: #" << length(needle_hits) << "\n";
      append(me.data_endPositions, getProperty(me.data_terminalStateMap, current_state)); // indices into TNeedle!
    }
  }


  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle, typename TWalker>
  inline bool _consumeChar(typename Pattern<TNeedle, AhoCorasickAmb>& me,
                           typename TWalker& walker, // a MasterVertex or Spawn
                           AminoAcid c,
                           typename Tag<FixedAASpec> fixedAASpec)
  {
    // if we cannot go down, but up is possible:
    while ((!goDown(me, walker, c)) &&
             goUp(me, walker)) // returns false when Spawn goes out of scope
    { /* .. follow suffix links upwards */ }
    if (atRoot(me, walker))
    {
      return false;
    }
    else // found a successor
    {
      addHits(me, walker);
      return true;
    }
  }

  // 
  // This is called by Spawns only!
  // Returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _consumeChar(typename Pattern<TNeedle, AhoCorasickAmb>& me,
                           Spawn<TNeedle>& current,
                           AminoAcid c)
  {
    // see if the Spawn can take another ambAA...
    if (isAmbiguous(c))
    {
      if (current.ambAA_seen >= me.max_ambAA) return false; // die -- do not even try to create sub-spawns
      else ++current.ambAA_seen; // increase ambAA count -- also for sub-Spawns which will follow from here...
      c = _createSpawns(me, current, c); // ... and leave fixed AA for master spawn
      return _consumeChar(me, current, c, Tag<FixedAASpec>());
    }
    return _consumeChar(me, current, c, Tag<FixedAASpec>());
  }

  // 
  // This is called by the MASTER thread only!
  // Returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _consumeChar(typename Pattern<TNeedle, AhoCorasickAmb>& me,
                           typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state,
                           AminoAcid c,
                           typename Tag<IsAmbAASpec> isAmbAASpec)
  {
    // if we reached the max# ambAA, we first need to go up, until the first ambAA goes out of scope
    if ((!me.ambAA_positions.empty()) &&
        (me.ambAA_positions.size() >= me.max_ambAA))
    {
      typedef Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
      KeyWordLengthType cum_depth_diff(0);
      Pattern<TNeedle, AhoCorasickAmb>::TVert old_state = current_state;
      while ((me.ambAA_positions.front() < cum_depth_diff) &&
             (goUp(me, current_state)))
      {
        cum_depth_diff += getProperty(me.data_nodeDepths, old_state) - getProperty(me.data_nodeDepths, current_state);
        old_state = current_state;
      }
      if (current_state == getRoot(me.data_graph))
      { // all the way to the top... reset all
        me.ambAA_positions.clear();
      }
      else
      { // update AA positions
        me.ambAA_positions.pop_front(); // first hit is invalid in any case
        typename Pattern<TNeedle, AhoCorasickAmb>::AmbAAPositionsIt it = me.ambAA_positions.begin();
        while (it != me.ambAA_positions.end())
        {
          if ((*it) < cum_depth_diff) // this position is out of scope
          { // remove and advance
            it = me.ambAA_positions.erase(it);
          }
          else
          { // update
            (*it) -= cum_depth_diff;
            ++it; // next AAPos
          }
        }
      }
    } // now, master can accept an ambAA again ...
    
    // then we add the current ambAA to the end
    me.ambAA_positions.push_back(getProperty(me.data_nodeDepths, current_state));

    c = _createSpawns(me, current_state, c); // ... and leave fixed AA for master thread
    return _consumeChar(me, current_state, c, Tag<FixedAASpec>());
  }

  // 
  // This is called by the master thread only!
  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _consumeChar(typename Pattern<TNeedle, AhoCorasickAmb>& me,
                           typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current,
                           AminoAcid c)
  {
    if (isAmbiguous(c)) 
      return _consumeChar(me, current, c, Tag<IsAmbAASpec>());
    else 
      return _consumeChar(me, current, c, Tag<FixedAASpec>());
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
  inline bool find(TFinder & finder, Pattern<TNeedle, AhoCorasickAmb> & me) {
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVert;

    if (empty(finder))
    {
      _patternInit(me);
      _finderSetNonEmpty(finder);
    }
    else
    {
      finder += me.data_needleLength; // restore last consumed position in haystack
      if (!empty(me.data_endPositions))
      { // Process left-over hits
        _reportHit(finder, me);
        return true;
      }
      ++finder; // advance to next position
    }

    while (!atEnd(finder))
    {
      const AminoAcid c = *finder;
      DEBUG_ONLY std::cout << "-- consuming " << c << " ---\n";
      // spawns; do them first, since we might add new spawns in main-thread & sub-spawns which are however settled at that point
      if (!me.spawns.empty()) {
        typename Pattern<TNeedle, AhoCorasickAmb>::SpawnIt it = me.spawns.begin();
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
      
      if (!empty(me.data_endPositions))
      {
        _reportHit(finder, me);
        return true;
      }

      ++finder;
    }
    return false;
  }

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H
