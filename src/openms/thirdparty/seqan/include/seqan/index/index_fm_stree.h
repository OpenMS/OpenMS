// ==========================================================================
//                 seqan - the library for sequence analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_STREE_H_
#define INDEX_FM_STREE_H_

namespace seqan {

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class VertexFmi
// ----------------------------------------------------------------------------

template <typename TAlphabet, typename TSize>
struct VertexFmi
{
    Pair<TSize> range;
    TSize       repLen;
    TAlphabet   lastChar;

    VertexFmi() :
        range(0, 0),
        repLen(0),
        lastChar(0)
    {}

    VertexFmi(MinimalCtor) :
        range(0, 0),
        repLen(0),
        lastChar(0)
    {}

    VertexFmi(Pair<TSize> newCurrentRange, TSize newRepLen, TAlphabet newChar) :
        range(newCurrentRange),
        repLen(newRepLen),
        lastChar(newChar)
    {}

    VertexFmi(VertexFmi const & other) :
        range(other.range),
        repLen(other.repLen),
        lastChar(other.lastChar)
    {}

    inline VertexFmi &
    operator = (VertexFmi const & _origin)
    {
        range = _origin.range;
        repLen = _origin.repLen;
        lastChar = _origin.lastChar;
        return *this;
    }
};

// ----------------------------------------------------------------------------
// Class HistoryStackFmi_
// ----------------------------------------------------------------------------

template <typename TAlphabet, typename TSize>
struct HistoryStackFmi_
{
    Pair<TSize> range;		// current SA interval of hits
    TAlphabet   lastChar;

    HistoryStackFmi_() {}

    template <typename TAlphabet_, typename TSize_>
    HistoryStackFmi_(TAlphabet_ const _lastChar, Pair<TSize_> const &_range):
        range(_range),
        lastChar(_lastChar)
    {}

    inline HistoryStackFmi_ const &
    operator=(HistoryStackFmi_ const & _origin)
    {
        range = _origin.range;
        lastChar = _origin.lastChar;
        return _origin;
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction VertexDescriptor                                        [Index]
// ----------------------------------------------------------------------------

template < typename TText, typename TOccSpec, typename TIndexSpec>
struct VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >
{
    typedef Index<TText,FMIndex<TOccSpec, TIndexSpec> >     TIndex;
    typedef typename Value<TIndex>::Type                    TAlphabet;
    typedef typename Size<TIndex>::Type                     TSize;

    typedef VertexFmi<TAlphabet, TSize> Type;
};

// ----------------------------------------------------------------------------
// Metafunction HistoryStackEntry_                                      [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TSpec, typename TIterSpec>
struct HistoryStackEntry_<Iter<Index<TText, FMIndex<TOccSpec, TSpec> >,
                               VSTree< TopDown< ParentLinks<TIterSpec> > > > >
{
    typedef HistoryStackFmi_<typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type,
                             typename Size<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type>     Type;
};

// ----------------------------------------------------------------------------
// Metafunction EdgeLabel                                            [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
struct EdgeLabel<Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > >
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _indexRequireTopDownIteration()                             [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec>
void _indexRequireTopDownIteration(Index<TText, FMIndex<TOccSpec, TIndexSpec> > & index)
{
    indexRequire(index, FibreSaLfTable());
}

// ----------------------------------------------------------------------------
// Function begin()                                                  [Iterator]
// ----------------------------------------------------------------------------

// ==========================================================================
///.Function.begin.param.object.type:Spec.FMIndex
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TIndexSpec> >, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TIndexSpec> > & index, TSpec const /*Tag*/)
{
	typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, TSpec>::Type TIter;

	TIter it(index);
	value(it).range.i1 = index.lfTable.prefixSumTable[0];

	return it;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TIndexSpec> > const, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TIndexSpec> > const & index, TSpec const /*Tag*/)
{
	typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> > const, TSpec>::Type TIter;

	TIter it(index);
	value(it).range.i1 = index.lfTable.prefixSumTable[0];

	return it;
}

// ----------------------------------------------------------------------------
// Function _isRoot()                                                [Iterator]
// ----------------------------------------------------------------------------

template <typename TAlphabet, typename TSize>
inline bool _isRoot(VertexFmi<TAlphabet, TSize> const & value)
{
    return _isSizeInval(value.range.i2);
}

// ----------------------------------------------------------------------------
// Function _isLeaf()                                                [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const & it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    return (value(it).range.i1 + 1 >= value(it).range.i2 &&
            value(it).range.i1 == _getSentinelPosition(container(it).lfTable.occTable));
}

template <typename TText, typename TSetSpec, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const & it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TSetSpec, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    return (value(it).range.i1 + 1 >= value(it).range.i2 &&
            sentinelPosition(getFibre(getFibre(container(it), FibreLfTable()), FibreOccTable()), value(it).range.i1));
}

// ----------------------------------------------------------------------------
// Function _getNodeByChar()                                         [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
inline bool _getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it,
                           typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type const & vDesc,
                           Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> & _range,
                           TChar c)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type                 TAlphabetSize;
    typedef typename Size<TIndex>::Type                         TSize;

    typedef typename Fibre<TIndex, FibreLfTable>::Type          TLfTable;
    typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

    TIndex const & _index = container(it);
    TPrefixSumTable const & pst = getFibre(getFibre(_index, FibreLfTable()), FibrePrefixSumTable());

    TAlphabetSize cPosition = getCharacterPosition(pst, c);

    if (_isRoot(vDesc))
    {
        _range.i1 = getPrefixSum(pst, cPosition);
        _range.i2 = getPrefixSum(pst, cPosition + 1);
    }
    else
    {
        TSize prefixSum = getPrefixSum(pst, cPosition);
        _range.i1 = prefixSum + countOccurrences(_index.lfTable.occTable, c, vDesc.range.i1 - 1);
        _range.i2 = prefixSum + countOccurrences(_index.lfTable.occTable, c, vDesc.range.i2 - 1);
    }

    return _range.i1 + 1 <= _range.i2;
}

// ----------------------------------------------------------------------------
// Function _goDownChar()                                            [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
inline bool _goDownChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, TChar c)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef Pair<typename Size<TIndex>::Type>                   TRange;

    TRange _range;

    if (_getNodeByChar(it, value(it), _range, c))
    {
        _historyPush(it);

        value(it).range = _range;
        value(it).lastChar = c;
        value(it).repLen++;

        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function _goDown()                                                [Iterator]
// ----------------------------------------------------------------------------

// TODO(esiragusa): Implement this.
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, FMIndex<TOccSpec,TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goDown(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TSpec, typename TIndexSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef typename Value<TIndex>::Type                    TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type             TAlphabetSize;

    if (isLeaf(it)) return false;

    for (TAlphabetSize c = 0; c < ValueSize<TAlphabet>::VALUE; ++c)
        if (_goDownChar(it, c)) return true;

    return false;
}

// ----------------------------------------------------------------------------
// Function _goDownString()                                          [Iterator]
// ----------------------------------------------------------------------------

//template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
//inline bool _goDownString(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
//                          TString const & string,
//                          TSize & lcp)
//{
//    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
//    typedef Pair<typename Size<TIndex>::Type>                   TRange;
//
//    typedef typename Iterator<TString const, Standard>::Type    TIterator;
//
//    _historyPush(it);
//
//    lcp = 0;
//
//    for (TIterator stringIt = begin(string, Standard()); stringIt != end(string, Standard()); ++stringIt)
//    {
//        TRange _range;
//
//        if (isLeaf(it) || !_getNodeByChar(it, value(it), _range, value(stringIt)))
//            return false;
//
//        value(it).range = _range;
//        value(it).lastChar = value(stringIt);
//        value(it).repLen++;
//        lcp++;
//    }
//
//    return true;
//}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool _goDownString(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                          TString const & string,
                          TSize & lcp)
{
    //typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    //typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename Iterator<TString const, Standard>::Type    TStringIter;

    lcp = 0;

    for (TStringIter stringIt = begin(string, Standard()); stringIt != end(string, Standard()); ++stringIt, ++lcp)
        if (isLeaf(it) || !_goDownChar(it, value(stringIt)))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function _goRight()                                               [Iterator]
// ----------------------------------------------------------------------------

// TODO(esiragusa): Implement this.
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goRight(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                     VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goRight(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goRight(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                     VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type                 TAlphabetSize;
    typedef Pair<typename Size<TIndex>::Type>                   TRange;

    typedef typename VertexDescriptor<TIndex>::Type             TVertexDescriptor;

    if (isRoot(it)) return false;

    TVertexDescriptor parentDesc = nodeUp(it);
    TRange _range;

    for (TAlphabetSize c = ordValue(value(it).lastChar) + 1; c < ValueSize<TAlphabet>::VALUE; ++c)
        if (_getNodeByChar(it, parentDesc, _range, c))
        {
            value(it).range = _range;
            value(it).lastChar = c;

            return true;
        }
    
    return false;
}

// ----------------------------------------------------------------------------
// Function _goUp()                                                  [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool _goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
	if (!isRoot(it))
	{
        value(it).range = it._parentDesc.range;
        value(it).lastChar = it._parentDesc.lastChar;
        --value(it).repLen;
        return true;
    }
    
    return false;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool _goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
	if (!isRoot(it))
	{
        value(it).range = back(it.history).range;
        value(it).lastChar = back(it.history).lastChar;
        --value(it).repLen;
        pop(it.history);
        return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------
// Function nodeUp()                                                 [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, class TSpec>
inline typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type
nodeUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > const & it)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef typename VertexDescriptor<TIndex>::Type         TVertexDescriptor;

    if (!empty(it.history))
        return TVertexDescriptor(back(it.history).range, value(it).repLen - 1, back(it.history).lastChar);
    else
        return value(it);
}

// ----------------------------------------------------------------------------
// Function _historyPush()                                           [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline void _historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > & it)
{
    it._parentDesc = value(it);
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline void
_historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    typedef Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;

    typename HistoryStackEntry_<TIter>::Type h;

    h.range = value(it).range;
    h.lastChar = value(it).lastChar;

    appendValue(it.history, h);
}

// ----------------------------------------------------------------------------
// Function repLength()                                              [Iterator]
// ----------------------------------------------------------------------------

template <typename TIndex, typename TAlphabet, typename TSize>
inline typename Size<TIndex>::Type
repLength(TIndex const &, VertexFmi<TAlphabet, TSize> const & vDesc)
{
	return vDesc.repLen;
}

// ----------------------------------------------------------------------------
// Function parentEdgeLabel()                                        [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline typename EdgeLabel<Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

// ----------------------------------------------------------------------------
// Function parentEdgeFirstChar()                                    [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline typename Value<Index<TText, FMIndex<TOccSpec, TIndexSpec > > >::Type
parentEdgeFirstChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

}
#endif  // INDEX_FM_STREE_H_
