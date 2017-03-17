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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

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
    KeyWordLengthType max_depth_decrease; // maximum loss in depths of traversed nodes (both while reporting hits and changing its own state)
    KeyWordLengthType ambAA_seen;         // number of ambAA's which the spawn has seen

    private:
	    Spawn();

    public:
	    Spawn(TVert init_state, KeyWordLengthType current_depth, KeyWordLengthType aaa_seen) :
		    current_state(init_state),
		    max_depth_decrease(current_depth),
        ambAA_seen(aaa_seen)
	    {}	
  };

  template <typename TNeedle>
  struct PatternHelperData
  {

    PatternHelperData()
    : data_endPositions(),
      data_keywordIndex(0),
      data_needleLength(0),
      data_lastState(0), // a bit of cheating, but we know that root==0
      spawns(),
      ambAA_positions()
    {}

    void reset()
    {
      clear(data_endPositions);
      data_keywordIndex = 0;
      data_needleLength = 0;
      data_lastState = 0; // a bit of cheating, but we know that root==0
      spawns.clear();
      ambAA_positions.clear();
    }

    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVert;
    typedef __uint8 KeyWordLengthType;

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
    typedef typename std::list<KeyWordLengthType>::const_iterator AmbAAPositionsCIt;
    AmbAAPositions ambAA_positions;    // indices of ambAA's relative to current path in trie; when going up, this list must be updated
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
    const KeyWordLengthType max_ambAA;            // default: 3
                                        
    // "constant" data, after construction of trie
	  Holder<TNeedle> data_host;                    // holds needles
    String<String<TSize> > data_terminalStateMap; // regular trie data -- plus: this gets augmented with all suffix traversals which are output nodes
    TGraph data_graph;                            // regular trie data
    String<TVert> data_supplyMap;                 // trie suffix links
    String<KeyWordLengthType> data_nodeDepths;    // depths of each graph node

    //____________________________________________________________________________
    Pattern() {}

    /**
      @brief Pattern Ctor with vector of needles, i.e. keywords/peptides.

      The vector @p nld must not be empty!

    */
    Pattern(TNeedle const& ndl, KeyWordLengthType max_AAA = 3)
      : nilVal(getNil<TVert>()),
        max_ambAA(max_AAA)
    {
      SEQAN_CHECKPOINT
      typedef typename Value<TNeedle>::Type TKeyword;
      typedef typename Value<TKeyword>::Type TAlphabet;
      
      for (TSize i = 0; i < length(ndl); ++i)
      {
        if (length(ndl[i]) > std::numeric_limits<KeyWordLengthType>::max())
        {
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, "Pattern<AhoCorasickAmb>(PeptideSet)", std::string("Input peptide to AhoCorasickAmb must NOT be longer than 255 chars!").c_str(), std::string(begin(ndl[i]), end(ndl[i])));
        }
        for (TSize j = 0; j < length(ndl[i]); ++j)
        {
          if (isAmbiguous(ndl[i][j])) // this check is important -- find() code below relies on no ambiguous chars being present!
          {
            throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, "Pattern<AhoCorasickAmb>(PeptideSet)", std::string("Input peptide to AhoCorasickAmb must NOT contain ambiguous amino acids ('B'/'Z'/'X')! Note: unknown AAs (e.g. 'U') will be converted to 'X' implicitly!").c_str(), std::string(begin(ndl[i]), end(ndl[i])));
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
    clear(me.data_terminalStateMap);
    _createAcTrie(me);

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

  //____________________________________________________________________________
  template <typename TNeedle>
  inline typename Size<TNeedle>::Type position(const PatternHelperData<TNeedle>& dh)
  {
    return dh.data_keywordIndex;
  }

  template <typename TFinder, typename TNeedle>
  inline void _reportHit(TFinder& finder, const Pattern<TNeedle, AhoCorasickAmb>& me, PatternHelperData<TNeedle>& dh) {
    size_t idx_endPosVec = length(dh.data_endPositions) - 1;
    dh.data_keywordIndex = dh.data_endPositions[idx_endPosVec];
    resize(dh.data_endPositions, idx_endPosVec); // pop last hit
    dh.data_needleLength = length(value(host(me), dh.data_keywordIndex)) - 1;
    finder -= dh.data_needleLength; // position finder at beginning of hit
    _setFinderLength(finder, dh.data_needleLength + 1);
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


  struct IsAmbAASpec {};
  struct FixedAASpec {};

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
    static const T ord_b = ordValue(AminoAcid('B'));
#if NDEBUG
#else
    assert(ord_b == 20); // make sure the table is ordered as we expect
    assert(ordValue(AminoAcid('Z')) == 21); // make sure the table is ordered as we expect
    assert(ordValue(AminoAcid('X')) == 22); // make sure the table is ordered as we expect
#endif
    idxFirst = jump[ordValue(c) - ord_b][1];
    idxLast = jump[ordValue(c) - ord_b][2];
    return AminoAcid(jump[ordValue(c) - ord_b][0]); // AA for main thread
  }

  // returns false if it reached the 'root', true otherwise
  // This will only be called by the master itself, and at this point, we can surely add another ambAA (goUp() was called before)!
  template <typename TNeedle>
  inline AminoAcid _createSpawns(const Pattern<TNeedle, AhoCorasickAmb>& me,
                                 PatternHelperData<TNeedle>& dh,
                                 AminoAcid c) // ALWAYS ambiguous!!!!
  {
	  DEBUG_ONLY std::cout << "found AAA: " << c << "\n";
    typedef typename Size<AminoAcid>::Type TSize;
    typedef typename Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    TSize idxFirst, idxLast;
    c = _getSpawnRange(idxFirst, idxLast, c);
    for (TSize i = idxFirst; i <= idxLast; ++i)
    {
      typedef typename Pattern<TNeedle, AhoCorasickAmb>::TVert TVert;
      TVert node_spawn = dh.data_lastState; // a potential spawn
      if (_consumeChar(me, dh, node_spawn, AminoAcid(i), Tag<FixedAASpec>()))
      {
        KeyWordLengthType node_depth = getProperty(me.data_nodeDepths, node_spawn); // depth at which the AA was consumed!
        // count how many ambAA positions from Master were skipped while going up
        const KeyWordLengthType removed_prefix_len = getProperty(me.data_nodeDepths, dh.data_lastState) - (node_depth - 1); // level ups before consuming (must be positive): 0..N
        DEBUG_ONLY std::cout << "  Spawn removed_prefix_len: " << int(removed_prefix_len) << "\n";
        KeyWordLengthType ambAA_seen(1); // spawn has seen 1 ambAA (this one), plus whatever was on the Master-path minus 'removed_prefix_len'
        typename PatternHelperData<TNeedle>::AmbAAPositions::reverse_iterator it = dh.ambAA_positions.rbegin();
        for (; it != dh.ambAA_positions.rend(); ++it)
        { // start at end (biggest position), down to smallest ambAA position. Stop when prefix was cut.
          if ((*it) > removed_prefix_len) ++ambAA_seen;
          else break; 
        }
        // push_front is paramount, since we might iterate over old spawns at this very moment
        // spawn gets the max_depth_decrease from what the master would have as first ambAA after moving
        if (it == dh.ambAA_positions.rbegin())
        { // all ambAA's were forgotten or there are none --> current depths counts (?)
          DEBUG_ONLY std::cout << "  No ambAA in current path. Using current depth after settling (since this is an ambAA)\n";
        } 
        else {
          node_depth = (*(--it)) - removed_prefix_len; // back to surviving ambAA; and its position minus the removed prefix
          DEBUG_ONLY std::cout << "  Updating ambAA of 1st Spawn: delta is now: " << int(node_depth) << "\n";
        }
        dh.spawns.push_front(Spawn<TNeedle>(node_spawn, node_depth, ambAA_seen));
        DEBUG_ONLY std::cout << "  1st Spawn from Master consuming '" << AminoAcid(i) << "' created at delta: " << int(dh.spawns.front().max_depth_decrease) << " AA-seen: " << int(dh.spawns.front().ambAA_seen) << "\n";
      }
    }
    return c;
  }
  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline AminoAcid _createSpawns(const Pattern<TNeedle, AhoCorasickAmb>& me,
                                 PatternHelperData<TNeedle>& dh,
                                 const Spawn<TNeedle>& spawn,
                                 AminoAcid c) // ALWAYS ambiguous!!!!
  {
	  DEBUG_ONLY std::cout << "trying to spawn on AAA: " << c << "\n";
    typedef typename Size<AminoAcid>::Type TSize;
    typedef typename Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    TSize idxFirst, idxLast;
    c = _getSpawnRange(idxFirst, idxLast, c);
    for (TSize i = idxFirst; i <= idxLast; ++i)
    {
      Spawn<TNeedle> spawn2 = spawn; // a potential spawn
      if (_consumeChar(me, dh, spawn2, AminoAcid(i), Tag<FixedAASpec>()))
      {
        // Spawn2 inherits the depths from its parent, since the master will also see this same AAA and spawn himself
        dh.spawns.push_front(spawn2);
        DEBUG_ONLY std::cout << "Spawn from Spawn '" << AminoAcid(i) << "' created at d: " << int(spawn2.max_depth_decrease) << " AA-seen: " << int(spawn2.ambAA_seen) << "\n";
      }
    }
    return c;
  }

  /// ###  go down  ####
  template<class TNeedle> inline bool goDown(const Pattern<TNeedle, AhoCorasickAmb>& me, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current, AminoAcid c)
  { 
    typename Pattern<TNeedle, AhoCorasickAmb>::TVert successor = getSuccessor(me.data_graph, current, c);
    if (successor == me.nilVal) return false;
    DEBUG_ONLY std::cout << "master/test matched '" << c << "'\n";
    current = successor;
    return true;
  }
  template<class TNeedle> inline bool goDown(const Pattern<TNeedle, AhoCorasickAmb>& me, Spawn<TNeedle>& spawn, AminoAcid c)
  {
    typename Pattern<TNeedle, AhoCorasickAmb>::TVert successor = getSuccessor(me.data_graph, spawn.current_state, c);
    if (successor == me.nilVal) return false;
    DEBUG_ONLY std::cout << "spawn matched '" << c << "' AA-seen: " << int(spawn.ambAA_seen) << "\n";
    spawn.current_state = successor;
    return true;
  }

  /// ###  go up  ####
  template<class TNeedle> inline bool goUp(const Pattern<TNeedle, AhoCorasickAmb>& me, Spawn<TNeedle>& spawn)
  {
    //if (atRoot(me, spawn)) return false; // cannot happen -- spawn would have died before
    const typename Pattern<TNeedle, AhoCorasickAmb>::TVert suffix_node = getProperty(me.data_supplyMap, spawn.current_state);
    // check if spawn is allowed to loose that many chars in front
    typedef typename Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    const KeyWordLengthType depthDiff = getProperty(me.data_nodeDepths, spawn.current_state) - getProperty(me.data_nodeDepths, suffix_node);
    if (spawn.max_depth_decrease <= depthDiff) {
      DEBUG_ONLY std::cout << "spawn died while going up (AAA out of scope)\n";
      spawn.current_state = getRoot(me.data_graph); // reset to root -- indicating failure!
      return false; // this spawn just threw away its reason of existance (i.e. the AAA). Die!
    }
    spawn.max_depth_decrease -= depthDiff;
    spawn.current_state = suffix_node; // no need to check for nilVal, since we cannot reach root (depths runs out before!)
    return true;
  }
  template<class TNeedle> inline bool goUp(const Pattern<TNeedle, AhoCorasickAmb>& me, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state)
  {
    if (atRoot(me, current_state)) return false;
    typename Pattern<TNeedle, AhoCorasickAmb>::TVert suffix_node = getProperty(me.data_supplyMap, current_state);
    if (suffix_node != me.nilVal) {
      current_state = suffix_node;
      return true;
    }
    return false;
  }

  template<class TNeedle> inline bool atRoot(const Pattern<TNeedle, AhoCorasickAmb>& me, Spawn<TNeedle>& spawn) {return spawn.current_state == getRoot(me.data_graph);}
  template<class TNeedle> inline bool atRoot(const Pattern<TNeedle, AhoCorasickAmb>& me, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state)
  {
    //DEBUG_ONLY std::cout << "master/test remains at root\n";
    return current_state == getRoot(me.data_graph);
  }

  template<class TNeedle> inline void addHits(const Pattern<TNeedle, AhoCorasickAmb>& me, PatternHelperData<TNeedle>& dh, Spawn<TNeedle>& spawn)
  {
    typedef typename Size<TNeedle>::Type TSize;
    String<TSize> needle_hits = getProperty(me.data_terminalStateMap, spawn.current_state);
    if (length(needle_hits) > 0)
    {
      int path_length = getProperty(me.data_nodeDepths, spawn.current_state); // == length of current path to spawn
      int unambiguous_suffix_length = path_length - spawn.max_depth_decrease; // == length of suffix peptide which does not contain AAA
      DEBUG_ONLY std::cout << "  spawn adding hits which are more than '" << unambiguous_suffix_length << "' chars long (thus contain the AAA).\n";

      // but only report those which contain the AAA
      for (int i = 0; i < (int)length(needle_hits); ++i)
      {
        int hit_length = (int)length(value(host(me), needle_hits[i]));
        if (hit_length <= unambiguous_suffix_length) break; // assumption: terminalStateMap is sorted by length of hits! ... uiuiui...
        DEBUG_ONLY std::cout << "  spawn hit: #" << i << "\n"; // value(value(host(me), needle_hits[i]))
        append(dh.data_endPositions, needle_hits[i]); // append hits which still contain the AAA
      }
    }
  }
  template<class TNeedle> inline void addHits(const Pattern<TNeedle, AhoCorasickAmb>& me, PatternHelperData<TNeedle>& dh, typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state)
  {
    typedef typename Size<TNeedle>::Type TSize;
    String<TSize> needle_hits = getProperty(me.data_terminalStateMap, current_state);
    if (length(needle_hits))
    {
      DEBUG_ONLY std::cout << "master/test hit total count: #" << length(needle_hits) << "\n";
      append(dh.data_endPositions, getProperty(me.data_terminalStateMap, current_state)); // indices into TNeedle!
    }
  }


  /**
     @brief Universal fixed char consumer. Works for Master and Spawns.

     ... using template specializations on 'TWalker' for inner function calls.

     @return False if it reached the 'root', true otherwise
   */
  template <typename TNeedle, typename TWalker>
  inline bool _consumeChar(const Pattern<TNeedle, AhoCorasickAmb>& me,
                           PatternHelperData<TNeedle>& dh,
                           TWalker& walker, // a MasterVertex or Spawn
                           AminoAcid c,
                           Tag<FixedAASpec> /*fixedAASpec*/)
  {
    DEBUG_ONLY std::cout << "consuming " << c<< " ";

    // if we cannot go down, but up is possible:
    while ((!goDown(me, walker, c)) &&
             goUp(me, walker)) // returns false when Spawn goes out of scope
    { /* .. follow suffix links upwards */ }
    if (atRoot(me, walker))
    {
      DEBUG_ONLY std::cout << "fail\n";
      return false;
    }
    else // found a successor
    {
      DEBUG_ONLY std::cout << "ok\n";
      addHits(me, dh, walker);
      return true;
    }
  }

  // 
  // This is called by Spawns only!
  // Returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _consumeChar(const Pattern<TNeedle, AhoCorasickAmb>& me,
                           PatternHelperData<TNeedle>& dh,
                           Spawn<TNeedle>& spawn,
                           AminoAcid c)
  {
    // see if the Spawn can take another ambAA...
    if (isAmbiguous(c))
    {
      if (spawn.ambAA_seen >= me.max_ambAA) return false; // we are at max and cannot consume more ambAA's; do not even try to create sub-spawns
      else ++spawn.ambAA_seen; // increase ambAA count -- also for sub-Spawns which will follow from here...
      c = _createSpawns(me, dh, spawn, c); // ... and leave fixed AA for master spawn
      return _consumeChar(me, dh, spawn, c, Tag<FixedAASpec>());
    }
    return _consumeChar(me, dh, spawn, c, Tag<FixedAASpec>());
  }


  // 
  // This is called by the master thread only!
  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _consumeChar(const Pattern<TNeedle, AhoCorasickAmb>& me,
                           PatternHelperData<TNeedle>& dh,
                           typename Pattern<TNeedle, AhoCorasickAmb>::TVert& current_state,
                           AminoAcid c)
  {
    typename Pattern<TNeedle, AhoCorasickAmb>::TVert old_state;
    typedef typename Pattern<TNeedle, AhoCorasickAmb>::KeyWordLengthType KeyWordLengthType;
    KeyWordLengthType cum_depth_diff(0);
    bool was_ambAA = isAmbiguous(c); // remember status, since we change 'c' and cannot query its ambiguity afterwards
    if (was_ambAA && (me.max_ambAA > 0)) // if max_ambAA==0, 'c' in {X,B,Z} will not be found in the trie using exact matching (since we forbid in needles), thus will topple back to root automatically, which is exactly what we want
    {
      // if we reached the max# ambAA, we first need to go up, until the first ambAA goes out of scope
      if ((!dh.ambAA_positions.empty()) &&
          (dh.ambAA_positions.size() >= me.max_ambAA))
      {
        old_state = current_state;
        // go up until first AA is out of scope
        while ((cum_depth_diff < dh.ambAA_positions.front()) &&
          (goUp(me, current_state)))
        {
          cum_depth_diff += getProperty(me.data_nodeDepths, old_state) - getProperty(me.data_nodeDepths, current_state);
          old_state = current_state;
        }
        if (current_state == getRoot(me.data_graph))
        { // all the way to the top... reset all
          dh.ambAA_positions.clear();
        }
        else
        { // update AA positions
          dh.ambAA_positions.pop_front(); // first hit is invalid in any case
          typename PatternHelperData<TNeedle>::AmbAAPositionsIt it = dh.ambAA_positions.begin();
          while (it != dh.ambAA_positions.end())
          {
            if ((*it) <= cum_depth_diff) // this position is out of scope
            { // remove and advance
              it = dh.ambAA_positions.erase(it);
            }
            else
            { // update
              (*it) -= cum_depth_diff;
              ++it; // next AAPos
            }
          }
        }
      } // now, master can accept an ambAA again ...

      c = _createSpawns(me, dh, c); // ... and leave fixed AA for master thread
    } // END ambiguous

    // 'c' is now the unambiguous char we consume
    
    // re-using some vars (old_state, cum_depth_diff) needed up there anyways
    old_state = current_state;
    bool was_consumed = _consumeChar(me, dh, current_state, c, Tag<FixedAASpec>());
    //
    // update previous ambAA char positions after we've moved within the tree
    //
    if (!dh.ambAA_positions.empty())
    {
      if (current_state == getRoot(me.data_graph))
      { // all the way to the top... reset all
        dh.ambAA_positions.clear();
      }
      else
      { // we are not at root, i.e. we went down exactly once!
        // count the upward steps, i.e. depth difference between before and after consuming the char: x*up() + 1*down()
        cum_depth_diff = getProperty(me.data_nodeDepths, old_state) - (getProperty(me.data_nodeDepths, current_state) - 1); // ranges from 0..N; cannot be negative!
        typename PatternHelperData<TNeedle>::AmbAAPositionsIt it = dh.ambAA_positions.begin();
        while (it != dh.ambAA_positions.end())
        {
          if ((*it) <= cum_depth_diff) // this position is out of scope
          { // remove and advance
            it = dh.ambAA_positions.erase(it);
          }
          else
          { // update
            (*it) -= cum_depth_diff;
            ++it; // next AAPos
          }
        }
      }
    } // END: ambAA_positions update

    if (was_ambAA && was_consumed) // was ambAA and we consumed it
    { // ... add position to the end
      dh.ambAA_positions.push_back(getProperty(me.data_nodeDepths, current_state));
    }

    return was_consumed;
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
  inline bool find(TFinder& finder, const Pattern<TNeedle, AhoCorasickAmb>& me, PatternHelperData<TNeedle>& dh)
  {
    if (empty(finder))
    {
      _finderSetNonEmpty(finder);
    }
    else
    {
      finder += dh.data_needleLength; // restore last consumed position in haystack
      if (!empty(dh.data_endPositions))
      { // Process left-over hits
        _reportHit(finder, me, dh);
        return true;
      }
      ++finder; // advance to next position
    }

    while (!atEnd(finder))
    {
      const AminoAcid c = *finder;
      DEBUG_ONLY std::cout << "\n\n-- consuming " << c << " ---\n";
      // spawns; do them first, since we might add new spawns in main-thread & sub-spawns which are however settled at that point
      if (!dh.spawns.empty()) {
        typename PatternHelperData<TNeedle>::SpawnIt it = dh.spawns.begin();
        DEBUG_ONLY std::cout << " --> Spawns (" << dh.spawns.size() << " alive):\n";
        while (it != dh.spawns.end())
        {
          if (!_consumeChar(me, dh, *it, c)) // might create new spawns
          { // spawn reached root --> kill it
            it = dh.spawns.erase(it); // remove and advance
            DEBUG_ONLY std::cout << " Killed spawn (" << dh.spawns.size() << " alive):\n";
          }
          else ++it; // next spawn
        }
      }
      // main thread
      DEBUG_ONLY std::cout << " --> Main (AA-seen: " << dh.ambAA_positions.size() << ", d: " << int(getProperty(me.data_nodeDepths, dh.data_lastState)) << ")\n";
      _consumeChar(me, dh, dh.data_lastState, c); // might create new spawns
      DEBUG_ONLY std::cout << "  <-- Main end (AA-seen: " << dh.ambAA_positions.size() << ", d: " << int(getProperty(me.data_nodeDepths, dh.data_lastState)) << ")\n";
      if (!empty(dh.data_endPositions))
      {
        _reportHit(finder, me, dh);
        return true;
      }

      ++finder;
    }
    return false;
  }

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H
