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
// $Authors: Chris Bielow, Sandro Andreotti $
// --------------------------------------------------------------------------


#ifndef SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H
#define SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H

#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>

namespace seqan
{

  //////////////////////////////////////////////////////////////////////////////
  // AhoCorasickAmbiguous
  //////////////////////////////////////////////////////////////////////////////

  /**
  .Spec.AhoCorasickAmbiguous:
  ..summary: Multiple exact string matching using Aho-Corasick.
  ..general:Class.Pattern
  ..cat:Searching
  ..signature:Pattern<TNeedle, AhoCorasickAmbiguous>
  ..param.TNeedle:The needle type, a string of keywords.
  ...type:Class.String
  ..remarks.text:The types of the keywords in the needle container and the haystack have to match.
  ..remarks.text:Matching positions do not come in order because we report beginning positions of matches.
  ..remarks.text:Likewise, if multiple keywords match at a given position no pre-specified order is guaranteed.
  ..example.text:The following example program searches for three needles ($queries$) in two haystack sequences ($db$) using the Aho-Corasick algorithm.
  ..example.file:demos/find/finder_aho_corasick.cpp
  ..example.remark:Note that you have to provide a String of Strings for the queries at the moment. This will be changed to a StringSet in the future.
  ..example.text:When executed, this program will create the following output.
  ..example.output:DB      POS     ENDPOS  TEXT
  0       0       4       MARD
  0       3       7       DPLY
  1       1       6       VGGGG
  1       6       9       AAA
  ..include:seqan/find.h
  */

  ///.Class.Pattern.param.TSpec.type:Spec.AhoCorasickAmbiguous

  struct AhoCorasickAmbiguous_;
  typedef Tag<AhoCorasickAmbiguous_> AhoCorasickAmbiguous;

  //////////////////////////////////////////////////////////////////////////////

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

    Holder<TNeedle> data_host; // holds needles
    String<String<TSize> > data_terminalStateMap; // regular trie data
    TGraph data_graph;                            // regular trie data
    String<TVertexDescriptor> data_supplyMap;     // trie suffix links

                                                  // To restore the automaton after a hit
    String<TSize> data_endPositions;	// All remaining keyword indices
    TSize data_keywordIndex;			// Current keyword that produced a hit
    TSize data_needleLength;			// Last length of needle to reposition finder
    TVertexDescriptor data_lastState;   // Last state in the trie

                                        //____________________________________________________________________________

    Pattern() {
    }

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
  inline void
    _createAcTrie(Pattern<TNeedle, AhoCorasickAmbiguous> & me)
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
    resizeVertexMap(me.data_graph, me.data_supplyMap);
    assignProperty(me.data_supplyMap, root, nilVal);

    // Bfs Traversal
    typedef typename Iterator<TGraph, BfsIterator>::Type TBfsIterator;
    TBfsIterator it(me.data_graph, root);
    //std::map<int, int> connectivity;
    goNext(it); // skip root
    for (;!atEnd(it);goNext(it)) {
      const GetValue<TBfsIterator>::Type itval = *it; // dereferencing *it will give an index into the property array!

                                                      //int edgecount = outDegree(me.data_graph, itval);
                                                      //++connectivity[edgecount];

      TVertexDescriptor parent = getProperty(parentMap, itval);
      TAlphabet sigma = getProperty(parentCharMap, itval);
      TVertexDescriptor down = getProperty(me.data_supplyMap, parent);
      while ((down != nilVal) &&
        (getSuccessor(me.data_graph, down, sigma) == nilVal))
      {
        down = getProperty(me.data_supplyMap, down);
      }
      if (down != nilVal) {
        assignProperty(me.data_supplyMap, itval, getSuccessor(me.data_graph, down, sigma));
        String<TPosition> endPositions = getProperty(me.data_terminalStateMap, getProperty(me.data_supplyMap, itval));
        if (!empty(endPositions)) {
          String<TPosition> endPositionsCurrent = getProperty(me.data_terminalStateMap, itval);
          typedef typename Iterator<String<TPosition>, Rooted >::Type TStringIterator;
          TStringIterator sit = begin(endPositions);
          for (;!atEnd(sit);goNext(sit)) {
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
  inline void
    setHost(Pattern<TNeedle, AhoCorasickAmbiguous> & me, TNeedle2 & needle)
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
    me.data_lastState = getRoot(me.data_graph);
  }


  //____________________________________________________________________________


  template <typename TNeedle>
  inline typename Size<TNeedle>::Type
    position(Pattern<TNeedle, AhoCorasickAmbiguous> & me)
  {
    return me.data_keywordIndex;
  }


  template <typename TFinder, typename TNeedle>
  inline void _consumeHit(TFinder & finder, Pattern<TNeedle, AhoCorasickAmbiguous> & me) {
    size_t idx_endPosVec = length(me.data_endPositions) - 1;
    me.data_keywordIndex = me.data_endPositions[idx_endPosVec];
    resize(me.data_endPositions, idx_endPosVec); // pop last hit
    me.data_needleLength = length(value(host(me), me.data_keywordIndex)) - 1;
    finder -= me.data_needleLength; // position finder at beginning of hit
    _setFinderLength(finder, me.data_needleLength + 1);
    _setFinderEnd(finder, position(finder) + length(finder)); // end of match within haystack
    return;
  }

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
        _consumeHit(finder, me);
        return true;
      }
      ++finder; // advance to next position
    }

    TVertexDescriptor successor;
    TVertexDescriptor& current = me.data_lastState;
    const TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    while (!atEnd(finder)) {
      // if we cannot go down, and up is possible:
      while (((successor = getSuccessor(me.data_graph, current, *finder)) == nilVal) &&
        (getProperty(me.data_supplyMap, current) != nilVal))
      { // .. follow suffix links upwards
        current = getProperty(me.data_supplyMap, current);
      }
      // found a successor:
      if (successor != nilVal) {
        current = successor;
        me.data_endPositions = getProperty(me.data_terminalStateMap, current); // indices into TNeedle!
        if (!empty(me.data_endPositions)) {
          _consumeHit(finder, me);
          return true;
        }
      }
      else { // no fitting successor and no way upwards --> we are at root
        current = getRoot(me.data_graph);
      }

      ++finder;
    }
    return false;
  }

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_AHOCORASICKAMBIGUOUS_H
