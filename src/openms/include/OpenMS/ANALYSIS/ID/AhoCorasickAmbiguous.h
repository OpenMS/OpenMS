// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Chris Bielow, Tobias Rausch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_AHOCORASICKAMBIGUOUS_H
#define OPENMS_ANALYSIS_ID_AHOCORASICKAMBIGUOUS_H

#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#ifdef NDEBUG
#define DEBUG_ONLY if (false)
#else
#define DEBUG_ONLY if (true)
#endif

// the SeqAn implementation comes first. To use the OpenMS interface, see below.

namespace seqan
{

  // we will use AAcid Tables from SeqAn 2.x here, since they are complete
  
  template <typename T = void>
  struct OMS_TranslateTableAAToChar_
  {
      static char const VALUE[27];
  };
  template <typename T>
  char const OMS_TranslateTableAAToChar_<T>::VALUE[27] =
  { // B = D|N,  Z = E|Q,  J = I|L
  // we re-order the aminoacids such that ambiguous AA's are consecutive (which saves effort during their enumeration)
    'A', // 00 Ala Alanine
    'Y', // 23 Tyr Tyrosine               !    1(Y)
    'C', // 02 Cys Cystine
    'D', // 03 Asp Aspartic Acid   // B
    'N', // 13 Asn Asparagine      // B   !    4(N)
    'F', // 05 Phe Phenylalanine       
    'G', // 06 Gly Glycine
    'H', // 07 His Histidine
    'I', // 08 Ile Isoleucine      // J
    'L', // 11 Leu Leucine         // J   !    9(L)
    'K', // 10 Lys Lysine                 !    10(K)
    'W', // 22 Trp Tryptophan			  !    11(W)
    'M', // 12 Met Methionine
    'O', // 14 Pyl Pyrrolysine            !    13(O)
    'P', // 15 Pro Proline                !    14(P)
    'E', // 04 Glu Glutamic Acid   // Z   !    15(E)
    'Q', // 16 Gln Glutamine       // Z
    'R', // 17 Arg Arginine
    'S', // 18 Ser Serine
    'T', // 19 Thr Threonine
    'U', // 20 Selenocystein
    'V', // 21 Val Valine
    'B', // 01 Aspartic Acid, Asparagine  $    22(B) // the AmbAA's need to be consecutive (B,J,Z,X)
    'J', // 09 Leucine, Isoleucine        $    23(J)
    'Z', // 24 Glutamic Acid, Glutamine   $
    'X', // 25 Unknown (25)
    '*'  // 26 Terminator
  };

  template <typename T = void>
  struct OMS_TranslateTableCharToAA_
  {
      static char const VALUE[256];
  };

  template <typename T>
  char const OMS_TranslateTableCharToAA_<T>::VALUE[256] =
  { // char --> amino acid (unsigned char with values from 0..26)
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //0
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //1
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  26,  25,  25,  25,  25,  25, //2
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //3

  //    ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,
      25,   0,  22,   2,   3,  15,   5,   6,   7,   8,  23,  10,   9,  12,   4,  13, //4

  //   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    ,
      14,  16,  17,  18,  19,  20,  21,  11,  25,   1,  24,  25,  25,  25,  25,  25, //5

  //    ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,
      25,   0,  22,   2,   3,  15,   5,   6,   7,   8,  23,  10,   9,  12,   4,  13, //6

  //   p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,    ,
      14,  16,  17,  18,  19,  20,  21,  11,  25,   1,  24,  25,  25,  25,  25,      //7

      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //8
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //9
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //10
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //11
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //12
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //13
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //14
      25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25  //15
  };

  template <typename T = void>
  struct OMS_TranslateTableByteToAA_
  {
    static char const VALUE[256];
  };

  template <typename T>
  char const OMS_TranslateTableByteToAA_<T>::VALUE[256] =
  {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, //0
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  25,  25,  25,  25,  25, //1
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //2
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //3
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //4
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //5
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //6
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //7
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //8
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //9
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //10
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //11
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //12
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //13
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //14
    25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25  //15
  };
  // --------------------------------------------------------------------------
  // Amino Acid
  // --------------------------------------------------------------------------
  struct AAcid_ {};
  typedef SimpleType<unsigned char, AAcid_> AAcid;

  template <> struct ValueSize<AAcid>
  {
    typedef uint8_t Type;
    static const Type VALUE = 27;
  };

  template <> struct BitsPerValue<AAcid>
  {
    typedef uint8_t Type;
    static const Type VALUE = 5;
  };

  inline AAcid unknownValueImpl(AAcid *)
  {
    static const AAcid _result = AAcid('X');
    return _result;
  }
  inline void assign(char & c_target, AAcid const & source)
  {
    c_target = OMS_TranslateTableAAToChar_<>::VALUE[source.value];
  }

  template <>
  struct CompareTypeImpl<AAcid, uint8_t>
  {
    typedef AAcid Type;
  };

  inline void assign(AAcid & target, uint8_t c_source)
  {
    target.value = OMS_TranslateTableByteToAA_<>::VALUE[c_source];
  }

  template <>
  struct CompareTypeImpl<AAcid, char>
  {
    typedef AAcid Type;
  };

  inline void assign(AAcid & target, char c_source)
  {
    target.value = OMS_TranslateTableCharToAA_<>::VALUE[(unsigned char)c_source];
  }

  typedef String<AAcid, Alloc<void> > AAString; // identical to seqan::Peptide, but without the misnomer (or how would you encode proteins? as seqan::Peptide?)
  typedef Iterator<AAString, Rooted>::Type AAStringIterator;

  //////////////////////////////////////////////////////////////////////////////

  struct FuzzyAC_;
  typedef Tag<FuzzyAC_> FuzzyAC;

  /// state of an AC spawn, operating on a trie
  template <typename TNeedle>
  struct Spawn
  {
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVert;
    typedef typename Pattern<TNeedle, FuzzyAC>::KeyWordLengthType KeyWordLengthType;
    TVert current_state;
    KeyWordLengthType max_depth_decrease; // maximum loss in depths of traversed nodes (both while reporting hits and changing its own state)
    KeyWordLengthType ambAA_seen;         // number of ambAA's which the spawn has seen

    public:
	    Spawn(TVert init_state, KeyWordLengthType current_depth, KeyWordLengthType aaa_seen) :
		    current_state(init_state),
		    max_depth_decrease(current_depth),
        ambAA_seen(aaa_seen)
	    {}	
      Spawn() = default;
      Spawn& operator=(const Spawn&) = default;
  };

  template <typename TNeedle>
  struct PatternAuxData
  {

    PatternAuxData()
    : hits_endPositions(),
      data_keywordIndex(0),
      data_needleLength(0),
      data_lastState(0), // a bit of cheating, but we know that root==0
      spawns()
    {}

    void reset()
    {
      clear(hits_endPositions);
      data_keywordIndex = 0;
      data_needleLength = 0;
      data_lastState = 0; // a bit of cheating, but we know that root==0
      spawns.clear();
    }

    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVert;
    typedef __uint8 KeyWordLengthType;

    // "working" set; changes with every hit
    String<TSize> hits_endPositions;	// All remaining keyword indices
    TSize data_keywordIndex;			// Current keyword that produced a hit
    TSize data_needleLength;			// Last length of needle to reposition finder
    TVert data_lastState;   // Last state of master instance in the trie
    typedef typename std::list<Spawn<TNeedle> > Spawns;
    typedef typename std::list<Spawn<TNeedle> >::iterator SpawnIt;
    typedef typename std::list<Spawn<TNeedle> >::const_iterator SpawnCIt;
    Spawns spawns;                      // spawn instances currently walking the tree
  };

  template <typename TNeedle>
  class Pattern<TNeedle, FuzzyAC>
  {
    //____________________________________________________________________________
  private:
    Pattern(Pattern const& other);
    Pattern const& operator=(Pattern const & other);
    //____________________________________________________________________________
  public:
    //typedef typename Size<TNeedle>::Type TSize; // defaults to uint64, but uint32 is enough...
    typedef uint32_t TSize; // max. number of peptides allowed (4 billion; checked in init())
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVert;
    typedef __uint8 KeyWordLengthType;

    // constant after C'tor; 
    const TVert nilVal;                           // NULL pointer for trie; e.g. returned when no successor is found
                                        
    // "constant" data, after construction of trie
    KeyWordLengthType max_ambAA;            // default: 3
    Holder<TNeedle> data_host;                    // holds needles, i.e. Peptides
    TGraph data_graph;                            // regular trie data
    String<String<TSize> > data_map_outputNodes; // regular trie data -- plus: this gets augmented with all suffix traversals which are output nodes
    String<KeyWordLengthType> data_node_depth;    // depths of each graph node

#ifndef NDEBUG
    String<TVert> parentMap;  ///< allows to find parent of each node
#endif

    //____________________________________________________________________________
    Pattern() 
      : nilVal(getNil<TVert>()),
        max_ambAA(0)
    {}

    /**
      @brief Pattern Ctor with vector of needles, i.e. keywords/peptides.

      The vector @p ndl must not be empty!

    */
    void init(TNeedle const& ndl, KeyWordLengthType max_AAA = 3)
    {
      SEQAN_CHECKPOINT
      max_ambAA = max_AAA;
      if (std::numeric_limits<TSize>::max() < length(ndl)) // check 4 billion peptide limit; use length(ndl) directly, since TSize might be too small
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, "Pattern<FuzzyAC>(PeptideSet)", std::string("Input contains more than 2^32 peptides. Cannot create trie.").c_str(), OpenMS::String(length(ndl)));
      }
      TSize ln = (TSize)length(ndl);
      for (TSize i = 0; i < ln; ++i)
      {
        if (length(ndl[i]) > std::numeric_limits<KeyWordLengthType>::max())
        {
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, "Pattern<FuzzyAC>(PeptideSet)", std::string("Input peptide to FuzzyAC must NOT be longer than 255 chars!").c_str(), std::string(begin(ndl[i]), end(ndl[i])));
        }
        TSize lni = (TSize)length(ndl[i]);
        for (TSize j = 0; j < lni; ++j)
        {
          if (isAmbiguous(ndl[i][j])) // this check is important -- find() code below relies on no ambiguous chars being present!
          {
            throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, "Pattern<FuzzyAC>(PeptideSet)", std::string("Input peptide to FuzzyAC must NOT contain ambiguous amino acids (B/J/Z/X)!").c_str(), std::string(begin(ndl[i]), end(ndl[i])));
          }
        }
      }
      setHost(*this, ndl);
    }

    ~Pattern()
    {
      SEQAN_CHECKPOINT
    }
    //____________________________________________________________________________
  };


  //////////////////////////////////////////////////////////////////////////////
  // Host Metafunctions
  //////////////////////////////////////////////////////////////////////////////

  template <typename TNeedle>
  struct Host< Pattern<TNeedle, FuzzyAC> >
  {
    typedef TNeedle Type;
  };

  template <typename TNeedle>
  struct Host< Pattern<TNeedle, FuzzyAC> const>
  {
    typedef TNeedle const Type;
  };


  //////////////////////////////////////////////////////////////////////////////
  // Functions
  //////////////////////////////////////////////////////////////////////////////
  /**
  @brief given an ambAA @p c, return a range of AA's (including @p idxLast) which need to be spawned.

  The range can be converted to AAcids via 'AAcid(T)'.
  */
  template<typename T>
  inline void _getSpawnRange(const AAcid c, T& idxFirst, T& idxLast)
  {
    // jump table:                  // start of spawns    // end of spawns (including)
    static const T jump[4][2] = { { ordValue(AAcid('D')), ordValue(AAcid('N')) },  // B = D,N
    { ordValue(AAcid('I')), ordValue(AAcid('L')) },  // J = I,L
    { ordValue(AAcid('E')), ordValue(AAcid('Q')) },  // Z = E,Q
    { 0                   , 21 } };// X = A..V
    static const T ord_b = ordValue(AAcid('B'));
#ifndef NDEBUG
    assert(ordValue(AAcid('N')) - ordValue(AAcid('D')) == 1); // must be neighbours
    assert(ordValue(AAcid('Q')) - ordValue(AAcid('E')) == 1); // must be neighbours
    assert(ordValue(AAcid('V')) == 21); // make sure the table is ordered as we expect
                                        // row of jump table:
    assert(ord_b == 22); // make sure the table is ordered as we expect
    assert(ordValue(AAcid('J')) == 23); // make sure the table is ordered as we expect
    assert(ordValue(AAcid('Z')) == 24); // make sure the table is ordered as we expect
    assert(ordValue(AAcid('X')) == 25); // make sure the table is ordered as we expect
#endif
    idxFirst = jump[ordValue(c) - ord_b][0];
    idxLast = jump[ordValue(c) - ord_b][1];
  }


  template <typename TNeedle>
  inline void _createAcTrie(Pattern<TNeedle, FuzzyAC> & me)
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
    createTrie(me.data_graph, me.data_map_outputNodes, host(me));

    // Create parent map
    String<TVert> parentMap;  ///< allows to find parent of each node
    String<TAlphabet> parentCharMap;      ///< allows to find character that led us to current node
    resizeVertexMap(me.data_graph, parentMap);
    resizeVertexMap(me.data_graph, parentCharMap);
    // init all nodes to Nil
    for (TPosition i = 0; i < length(parentMap); ++i) {
      assignProperty(parentMap, i, nilVal);
    }
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itEd(me.data_graph);
    //int count_prop = 0;
    for (; !atEnd(itEd); goNext(itEd)) {
      //             property,  vertex            , value
      assignProperty(parentMap, targetVertex(itEd), sourceVertex(itEd));
      assignProperty(parentCharMap, targetVertex(itEd), label(itEd));
      //++count_prop;
    }
    //std::cout << "\nassigned " << count_prop << " props to " << length(parentMap) << " values\n";
#ifndef NDEBUG    
    me.parentMap = parentMap;  ///< allows to find parent of each node
#endif

    // Build AC
    TVert root = getRoot(me.data_graph);
    // properties....
    String<TVert> data_map_failurelink;                 // trie suffix links (temporary only)
    resizeVertexMap(me.data_graph, data_map_failurelink);  
    assignProperty(data_map_failurelink, root, nilVal);

    resizeVertexMap(me.data_graph, me.data_node_depth);  // node depths
    assignProperty(me.data_node_depth, root, 0);

    // Bfs Traversal
    typedef typename Iterator<TGraph, BfsIterator>::Type TBfsIterator;
    TBfsIterator it(me.data_graph, root);
    typedef typename Size<AAcid>::Type TSize;
    TSize idxAAFirst, idxAALast; // range of unambiguous AAcids: AAcid(idx)
    _getSpawnRange('X', idxAAFirst, idxAALast);
    // create nextMove function for root (point to itself)
    for (TSize idx = idxAAFirst; idx <= idxAALast; ++idx)
    {
      if (getSuccessor(me.data_graph, root, AAcid(idx)) == nilVal) addEdge(me.data_graph, root, root, AAcid(idx));
    }

    goNext(it); // skip root

    // connectivity statistics
    //std::map<int, int> connectivity;

    for (; !atEnd(it); goNext(it))
    {
      const typename GetValue<TBfsIterator>::Type itval = *it; // dereferencing *it will give an index into the property array!

      //++connectivity[outDegree(me.data_graph, itval)];

      // set depth of current node using: depths(parent) + 1
      TVert parent = getProperty(parentMap, itval);
      assignProperty(me.data_node_depth, itval, getProperty(me.data_node_depth, parent) + 1);


      ///
      /// create failure function (suffix links) and output function
      ///
      // sigma: edge label
      TAlphabet sigma = getProperty(parentCharMap, itval);
      // take suffix link of parent and try to go down with sigma
      TVert down = getProperty(data_map_failurelink, parent);
      while ((down != nilVal) &&
        (getSuccessor(me.data_graph, down, sigma) == nilVal))
      {
        down = getProperty(data_map_failurelink, down);
      }
      if (down != nilVal)
      { // we found an edge to follow down
        assignProperty(data_map_failurelink, itval, getSuccessor(me.data_graph, down, sigma));
        // output function
        String<TPosition> endPositions = getProperty(me.data_map_outputNodes, getProperty(data_map_failurelink, itval));
        if (!empty(endPositions))
        {
          // get current end positions (full path) ...
          String<TPosition> endPositionsCurrent = getProperty(me.data_map_outputNodes, itval);
          typedef typename Iterator<String<TPosition>, Rooted >::Type TStringIterator;
          // .. and append all patterns which are a suffix
          TStringIterator sit = begin(endPositions);
          for (;!atEnd(sit); goNext(sit))
          {
            appendValue(endPositionsCurrent, *sit);
          }
          assignProperty(me.data_map_outputNodes, itval, endPositionsCurrent);
        }
      }
      else { // no suffix exists: point suffix link of current node to root
        assignProperty(data_map_failurelink, itval, root);
      }
      
      // create nextMove function
      for (TSize idx = idxAAFirst; idx <= idxAALast; ++idx)
      {
        if (getSuccessor(me.data_graph, itval, AAcid(idx)) == nilVal)
        { // no child:
          const TVert& target = getSuccessor(me.data_graph, getProperty(data_map_failurelink, itval) , AAcid(idx));
          addEdge(me.data_graph, itval, target, AAcid(idx));
        }
        
      }
    }

    //sw.stop();
    //std::cout << OpenMS::String(" Trie build time: ") + sw.getClockTime() + " s (wall), " + sw.getCPUTime() + " s (CPU)).\n\n";

    /*for (std::map<int,int>::const_iterator it = connectivity.begin(); it != connectivity.end(); ++it)
    {
      std::cout << "   conn: " << it->first << " : " << it->second << "x\n";
    }*/

  }

  template <typename TNeedle, typename TNeedle2>
  void setHost(Pattern<TNeedle, FuzzyAC> & me, TNeedle2 const & needle)
  {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT(empty(needle));
    setValue(me.data_host, needle);
    clear(me.data_graph);
    clear(me.data_map_outputNodes);
    _createAcTrie(me);

    //fstream strm;
    //strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
    //String<String<char> > nodeMap;
    //_createTrieNodeAttributes(me.data_graph, me.data_map_outputNodes, nodeMap);
    //String<String<char> > edgeMap;
    //_createEdgeAttributes(me.data_graph,edgeMap);
    //write(strm,me.data_graph,nodeMap,edgeMap,DotDrawing());
    //strm.close();
    // Supply links
    //for(unsigned int i=0;i<length(me.data_map_failurelink);++i) {
    //	std::cout << i << "->" << getProperty(me.data_map_failurelink,i) << ::std::endl;
    //}
  }
  template <typename TNeedle, typename TNeedle2>
  inline void setHost(Pattern<TNeedle, FuzzyAC> & me, TNeedle2 & needle)
  {
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
  }
  //____________________________________________________________________________

  //____________________________________________________________________________
  template <typename TNeedle>
  inline typename Size<TNeedle>::Type position(const PatternAuxData<TNeedle>& dh)
  {
    return dh.data_keywordIndex;
  }

  template <typename TFinder, typename TNeedle>
  inline void _reportHit(TFinder& finder, const Pattern<TNeedle, FuzzyAC>& me, PatternAuxData<TNeedle>& dh) {
    size_t idx_endPosVec = length(dh.hits_endPositions) - 1;
    dh.data_keywordIndex = dh.hits_endPositions[idx_endPosVec];
    resize(dh.hits_endPositions, idx_endPosVec); // pop last hit
    dh.data_needleLength = length(value(host(me), dh.data_keywordIndex)) - 1;
    finder -= dh.data_needleLength; // position finder at beginning of hit
    _setFinderLength(finder, dh.data_needleLength + 1);
    _setFinderEnd(finder, position(finder) + length(finder)); // end of match within haystack
    return;
  }

  inline bool isAmbiguous(AAcid c)
  { // relies on the fact that ambiguous AA's are occuring in a block
    static const typename ValueSize<AAcid>::Type vB = ordValue(AAcid('B')); // D,N
    static const typename ValueSize<AAcid>::Type vX = ordValue(AAcid('X')); // all
    return (vB <= ordValue(c) && ordValue(c) <= vX);
  }

  inline bool isAmbiguous(AAString s)
  {
    for (AAStringIterator it = begin(s); it != end(s); ++it)
    {
      if (isAmbiguous(*it)) return true;
    }
    return false;
  }



  // This will only be called by the master itself (before moving to root); so far, no AAA was consumed (by master)
  template <typename TNeedle>
  inline void _createPrimarySpawns(const Pattern<TNeedle, FuzzyAC>& me,
                                   PatternAuxData<TNeedle>& dh,
                                   const AAcid c) // ALWAYS ambiguous!!!!
  {
    assert(isAmbiguous(c));

    DEBUG_ONLY std::cout << "found AAA: " << c << "\n";
    typedef typename Size<AAcid>::Type TSize;
    typedef typename Pattern<TNeedle, FuzzyAC>::TVert TVert;
    TSize idxFirst, idxLast;
    _getSpawnRange(c, idxFirst, idxLast);
    for (; idxFirst <= idxLast; ++idxFirst)
    {
      TVert node_spawn = dh.data_lastState; // last state of master
      if (_consumeChar(me, dh, node_spawn, AAcid(idxFirst))) // call this using master's _consumeChar(), since it might pass through root (which is allowed), but should not die.
      { // spawn from current position; push front to flag as 'processed' for the current input char
        // depths is 'current_depth - 1' (must be computed here!); ambAA-count: fixed to 1 (first AAA, since spawned from master)
        dh.spawns.push_front(Spawn<TNeedle>(node_spawn, getProperty(me.data_node_depth, node_spawn) - 1, 1));
        DEBUG_ONLY std::cout << "  Init Spawn from Master consuming '" << AAcid(idxFirst) << "\n";
      }
    }
  }

  // called by spawns only
  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _createSecondarySpawns(const Pattern<TNeedle, FuzzyAC>& me,
                                     PatternAuxData<TNeedle>& dh,
                                     Spawn<TNeedle>& spawn,
                                     const AAcid c) // ALWAYS ambiguous!!!!
  {
    assert(isAmbiguous(c));

    DEBUG_ONLY std::cout << "\n\ntrying to spawn from spawn on AAA: " << c << " with path: " << getPath(me, spawn.current_state) << "\n";
    typedef typename Size<AAcid>::Type TSize;
    TSize idxFirst, idxLast;
    _getSpawnRange(c, idxFirst, idxLast);
    for (TSize idx = idxFirst + 1; idx <= idxLast; ++idx) // first iteration is left for the parent
    {
      Spawn<TNeedle> spawn2 = spawn; // a potential spawn
      //std::cout << "spawn aa is " << AAcid(idx) << "\n";
      if (_consumeChar(me, dh, spawn2, AAcid(idx)))
      {
        // Spawn2 inherits the depths from its parent
        dh.spawns.push_front(spawn2);
        DEBUG_ONLY std::cout << "  Spawn from Spawn '" << getPath(me, spawn2.current_state) << "' created at d: " << int(spawn2.max_depth_decrease) << " AA-seen: " << int(spawn2.ambAA_seen) << "\n";
      }
    }
    bool r = _consumeChar(me, dh, spawn, AAcid(idxFirst));
    if (r)
    {
      DEBUG_ONLY std::cout << "  Spawn from Spawn '" << getPath(me, spawn.current_state) << "' created at d: " << int(spawn.max_depth_decrease) << " AA-seen: " << int(spawn.ambAA_seen) << "\n";
    }
    return r;
  }

#ifdef NDEBUG  
  template<class TNeedle> inline std::string getPath(const Pattern<TNeedle, FuzzyAC>& /*me*/, typename Pattern<TNeedle, FuzzyAC>::TVert /*current_state*/)
  {
    return "";
  }
#else
  /// for debug only
  template<class TNeedle> inline std::string getPath(const Pattern<TNeedle, FuzzyAC>& me, typename Pattern<TNeedle, FuzzyAC>::TVert current_state)
  {
    if (getRoot(me.data_graph) == current_state) return "";

    typename Pattern<TNeedle, FuzzyAC>::TVert suffix_node = getProperty(me.parentMap, current_state);
    typename Iterator<typename Pattern<TNeedle, FuzzyAC>::TGraph const, OutEdgeIterator>::Type it(me.data_graph, suffix_node);
    while(targetVertex(it) != current_state) ++it;
    char c = (label(it));
    return getPath(me, suffix_node) + c;
  }
#endif

  template<class TNeedle> inline void addHits(const Pattern<TNeedle, FuzzyAC>& me, PatternAuxData<TNeedle>& dh, Spawn<TNeedle>& spawn)
  {
  // TODO check if at root and return immediately?
    typedef typename Pattern<TNeedle, FuzzyAC>::TSize TSize;
    const String<TSize>& needle_hits = getProperty(me.data_map_outputNodes, spawn.current_state);
    //DEBUG_ONLY std::cout << "spawn at path: " << getPath(me, spawn.current_state) << "\n";
    if (length(needle_hits))
    {
      int path_length = getProperty(me.data_node_depth, spawn.current_state); // == length of current path to spawn
      int unambiguous_suffix_length = path_length - spawn.max_depth_decrease; // == length of suffix peptide which does not contain AAA
      DEBUG_ONLY std::cout << "  spawn adding hits which are at least " << unambiguous_suffix_length << " chars long (thus contain the AAA).\n";

      // but only report those which contain the AAA
      for (auto it = begin(needle_hits); it != end(needle_hits); ++it)
      {
        int hit_length = (int)length(value(host(me), *it));
        if (hit_length < unambiguous_suffix_length) break; // assumption: terminalStateMap is sorted by length of hits! ... uiuiui...
        DEBUG_ONLY std::cout << "  add spawn hit: needle #" << *it << " as " << (value(host(me), *it)) << "\n";
        append(dh.hits_endPositions, *it); // append hits which still contain the AAA
      }
    }
  }
  template<class TNeedle> inline void addHits(const Pattern<TNeedle, FuzzyAC>& me, PatternAuxData<TNeedle>& dh, typename Pattern<TNeedle, FuzzyAC>::TVert& current_state)
  {
    typedef typename Pattern<TNeedle, FuzzyAC>::TSize TSize;
    //DEBUG_ONLY std::cout << "master at path: " << getPath(me, current_state) << "\n";
    const String<TSize>& needle_hits = getProperty(me.data_map_outputNodes, current_state);
    if (length(needle_hits))
    {
      DEBUG_ONLY std::cout << "master's new hits: total " << length(needle_hits) << " hits\n";
      append(dh.hits_endPositions, getProperty(me.data_map_outputNodes, current_state)); // indices into TNeedle!
    }
  }


  /**
     @brief Universal fixed char consumer. Works for Master and Spawns.

     ... using template specializations on 'TWalker' for inner function calls.
     Adds hits to @p dh after settling.

     @return False if it reached the 'root', true otherwise
   */
  template <typename TNeedle>
  inline bool _consumeChar(const Pattern<TNeedle, FuzzyAC>& me,
                           PatternAuxData<TNeedle>& dh,
                           typename Pattern<TNeedle, FuzzyAC>::TVert& current,
                           const AAcid c)
  {
    //DEBUG_ONLY std::cout << "consuming real char " << c << " ";
    current = getSuccessor(me.data_graph, current, c);
    assert(current != me.nilVal);
    if (current != getRoot(me.data_graph)) {
      addHits(me, dh, current);
      return true;
    }
    return false;
  }

  template <typename TNeedle>
  inline bool _consumeChar(const Pattern<TNeedle, FuzzyAC>& me,
    PatternAuxData<TNeedle>& dh,
    Spawn<TNeedle>& spawn,
    const AAcid c)
  {
    //DEBUG_ONLY std::cout << "consuming real char " << c << " ";
    typename Pattern<TNeedle, FuzzyAC>::TVert successor = getSuccessor(me.data_graph, spawn.current_state, c);
    assert(successor != me.nilVal);
    // check if prefix was lost
    if (getProperty(me.data_node_depth, spawn.current_state) >= getProperty(me.data_node_depth, successor))
    { // went at least one level up (and maybe one down again, hence equality)
      typedef typename Pattern<TNeedle, FuzzyAC>::KeyWordLengthType KeyWordLengthType;
      const KeyWordLengthType up_count = 1 + getProperty(me.data_node_depth, spawn.current_state) - getProperty(me.data_node_depth, successor);
      // the +1 is not valid if we end up at root (but there, the spawn is dead anyway)
      if (up_count > spawn.max_depth_decrease)
      {
        //DEBUG_ONLY std::cout << "spawn died while going up (AAA out of scope)\n";
        //spawn.current_state = getRoot(me.data_graph); // reset to root -- not required
        return false; // this spawn just threw away its reason of existance (i.e. the AAA). Die!
      }
      spawn.max_depth_decrease -= up_count;
    }
    spawn.current_state = successor;
    if (spawn.current_state != getRoot(me.data_graph)) {
      addHits(me, dh, spawn); // use spawn version for length checking!
      return true;
    }
    return false;
  }

  // 
  // This is called by Spawns only!
  // Returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline bool _spawnConsumeChar(const Pattern<TNeedle, FuzzyAC>& me,
                           PatternAuxData<TNeedle>& dh,
                           Spawn<TNeedle>& spawn,
                           const AAcid c)
  {
    // see if the Spawn can take another ambAA...
    if (isAmbiguous(c))
    {
      if (spawn.ambAA_seen >= me.max_ambAA) return false; // we are at max and cannot consume more ambAA's; do not even try to create sub-spawns
      else ++spawn.ambAA_seen; // increase ambAA count -- also for sub-Spawns which will follow from here...
      return _createSecondarySpawns(me, dh, spawn, c); // child spawns are preprended to spawn list; returns if main spawn survived
    }
    // UNambiguous
    return _consumeChar(me, dh, spawn, c);
  }

  // This is called by the master thread only!
  // returns false if it reached the 'root', true otherwise
  template <typename TNeedle>
  inline void _masterConsumeChar(const Pattern<TNeedle, FuzzyAC>& me,
                           PatternAuxData<TNeedle>& dh,
                           const AAcid c)
  {
    if (isAmbiguous(c))
    {
      if (me.max_ambAA > 0)
      {
        _createPrimarySpawns(me, dh, c);
      }
      else
      { // do not spawn anything
        DEBUG_ONLY std::cout << " --> Main found AAA, but none allowed. Resetting to root. No spawns created.\n";
      }
      dh.data_lastState = getRoot(me.data_graph); // reset Master to root
      return;
    } // end: ambiguous

    // 'c' is UN-ambiguous
    _consumeChar(me, dh, dh.data_lastState, c); // adds hits as well
  }


  template <typename TFinder, typename TNeedle>
  inline bool find(TFinder& finder, const Pattern<TNeedle, FuzzyAC>& me, PatternAuxData<TNeedle>& dh)
  {
    if (empty(finder))
    {
      _finderSetNonEmpty(finder);
    }
    else
    {
      finder += dh.data_needleLength; // restore last consumed position in haystack
      if (!empty(dh.hits_endPositions))
      { // Process left-over hits
        _reportHit(finder, me, dh);
        return true;
      }
      // nothing to report at this point
      ++finder; // advance to next position
    }

    while (!atEnd(finder))
    {
      const AAcid c = *finder;
      DEBUG_ONLY std::cout << "\n\n-- consuming " << c << " ---\n";
      // spawns; do them first, since we might add new (but settled) spawns in main-thread & sub-spawns
      if (!dh.spawns.empty())
      {
        typename PatternAuxData<TNeedle>::SpawnIt it = dh.spawns.begin();
        //DEBUG_ONLY std::cout << " --> Spawns (" << dh.spawns.size() << " alive):\n";
        while (it != dh.spawns.end())
        {
          if (!_spawnConsumeChar(me, dh, *it, c)) // might create new spawns
          { // spawn reached root --> kill it
            it = dh.spawns.erase(it); // remove and advance
            //DEBUG_ONLY std::cout << " Killed spawn (" << dh.spawns.size() << " alive):\n";
          }
          else
          { // next spawn
            ++it;
          }
        }
      }
      // main thread
      DEBUG_ONLY std::cout << " --> Main; d: " << int(getProperty(me.data_node_depth, dh.data_lastState)) << ")\n";
      _masterConsumeChar(me, dh, c); // might create new spawns
      DEBUG_ONLY std::cout << "  <-- Main end; d: " << int(getProperty(me.data_node_depth, dh.data_lastState)) << ")\n";

      // print current states
      DEBUG_ONLY std::cout << " --> POST: Main state: " << getPath(me, dh.data_lastState) << "\n";
      for (typename PatternAuxData<TNeedle>::SpawnIt it = dh.spawns.begin(); it != dh.spawns.end(); ++it)
      {
        DEBUG_ONLY std::cout << " --> POST: Spawn state: " << getPath(me, it->current_state) << "\n";
      }

      if (!empty(dh.hits_endPositions))
      {
        _reportHit(finder, me, dh);
        return true;
      }

      ++finder;
    }
    return false;
  }

}// namespace SEQAN_NAMESPACE_MAIN

namespace OpenMS
{

  //////////////////////////////////////////////////////////////////////////////
  // FuzzyAC
  //////////////////////////////////////////////////////////////////////////////

  /**
  @brief Extended Aho-Corasick algorithm capable of matching ambiguous amino acids in the pattern (i.e. proteins).

  ...
  Features:
  + blazingly fast
  + low memory usage
  + number of allowed ambAA's can be capped by user (default 3).


  This implementation is based on the original AC in SeqAn.

  */

  class AhoCorasickAmbiguous
  {
  public:
    typedef typename ::seqan::StringSet<::seqan::AAString> PeptideDB;
    typedef typename ::seqan::Pattern<PeptideDB, ::seqan::FuzzyAC> FuzzyACPattern;

    /**
      @brief Construct a trie from a set of peptide sequences (which are to be found in a protein).

      Peptides must not contain ambiguous characters (exception thrown otherwise) or unknown characters (such as J or U).
      Ambiguous characters are only allowed in protein sequences.
      
      Usage:
      Build the pattern only once and use it multiple times when running findNext().

      @param pep_db Set of peptides
      @param aaa_max Maximum of allowed ambiguous characters in the matching protein sequence
      @param pattern The pattern to be created
      @throws Exception::InvalidValue if a peptide contains an unknown (U,J,...) or ambiguous character
    */
    static void initPattern(const PeptideDB& pep_db, const int aaa_max, FuzzyACPattern& pattern)
    {
      pattern.init(pep_db, KeyWordLengthType(aaa_max));
    }

    /**
      @brief Default Ctor; call setProtein() before using findNext().

    */
    AhoCorasickAmbiguous()
      : finder_(),
      protein_(),
      dh_()
    {
    }
    /**
      @brief Prepare to start searching for hits in a new protein sequence.

      This only sets the sequence. No computation is performed. Use findNext() to enumerate the hits.

      @param protein_sequence Sequence (ambiguous characters allowed)
    */
    AhoCorasickAmbiguous(const String& protein_sequence)
      : finder_(),
        protein_(),
        dh_()
    {
      setProtein(protein_sequence);
    }

    /**
      @brief Reset to new protein sequence. All previous data is forgotten.
    */
    void setProtein(const String& protein_sequence)
    {
      protein_ = protein_sequence.c_str(); // we need an internal copy, since finder_ only keeps a pointer
      finder_ = ::seqan::Finder<seqan::AAString>(protein_);
      dh_.reset();
    }

    /**
      @brief Enumerate hits.

      @param pattern The pattern (i.e. trie) created with initPattern().
      @return False if end of protein is reached. True if a hit is found.
    */
    bool findNext(const FuzzyACPattern& pattern)
    {
      return ::seqan::find(finder_, pattern, dh_);
    }

    /**
      @brief Get index of hit into peptide database of the pattern.
      
      Only valid if findNext() returned true before.
    */
    Size getHitDBIndex()
    {
      return ::seqan::position(dh_);
    }

    /**
      @brief Offset into protein sequence where hit was found.

      Only valid if findNext() returned true before.
    */
    Int getHitProteinPosition()
    {
      return (Int)position(finder_);
    }

  private:
    typedef typename FuzzyACPattern::KeyWordLengthType KeyWordLengthType;

    // member
    ::seqan::Finder<seqan::AAString> finder_; ///< locate the next peptide hit in protein
    ::seqan::AAString protein_;               ///< the protein sequence - we need to store it since the finder only keeps a pointer to protein when constructed
    ::seqan::PatternAuxData<PeptideDB> dh_;   ///< auxilliary data to hold a state after searching
  }; // class FuzzyAC

} // namespace OpenMS

#endif //#ifndef OPENMS_ANALYSIS_ID_AHOCORASICKAMBIGUOUS_H
