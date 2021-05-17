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

#ifndef INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H
#define INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================
template <typename TOccTable, typename TPrefixSumTable>
struct LfTable;

template <typename TValue>
struct WaveletTree;

template<typename TValue> 
class RankDictionary;

template<typename TRankDictionarySpec, typename TSpec> 
class SentinelRankDictionary;

struct FibreTreeStructure_;
typedef Tag<FibreTreeStructure_> const FibreTreeStructure;

template <typename TChar, typename TSpec = void>
class RightArrayBinaryTree;

// ==========================================================================
// Tags 
// ==========================================================================
//
/**
.Tag.RightArrayBinaryTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.RightArrayBinaryTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a RightArrayBinaryTree.
..cat:RightArrayBinaryTree
..tag.FibreTreeStructureEncoding:The string encoding the wavelet tree structure.
..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

///.Metafunction.Fibre.param.TContainer.type:Class.RightArrayBinaryTree
///.Metafunction.Fibre.param.TSpec.type:Tag.RightArrayBinaryTree Fibres
struct FibreTreeStructureEncoding_;
typedef Tag<FibreTreeStructureEncoding_> const FibreTreeStructureEncoding;

// ==========================================================================
// Metafunctions
// ==========================================================================
template <typename TChar, typename TSpec>
struct Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>
{
    // TODO (singer): Description why we need + 2
    typedef typename BitVector_<Log2<ValueSize<TChar>::VALUE + 2>::VALUE>::Type TPos;
    typedef String<Pair<TChar, TPos> > Type;
};

// template <typename TChar, typename TSpec>
// struct Fibre<RightArrayBinaryTree<TChar, TSpec> const, FibreTreeStructureEncoding>
// {
//     typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type const Type;
// };

// ==========================================================================
template <typename TChar, typename TSpec>
struct Reference<RightArrayBinaryTree<TChar, TSpec> >
{
    typedef typename Value<RightArrayBinaryTree<TChar, TSpec> >::Type & Type;
};

template <typename TChar, typename TSpec>
struct Reference<const RightArrayBinaryTree<TChar, TSpec> >
{
    typedef typename Value<RightArrayBinaryTree<TChar, TSpec> >::Type const & Type;
};

// ==========================================================================
template <typename TChar, typename TSpec>
struct Value<RightArrayBinaryTree<TChar, TSpec> >
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type TWaveletTreeVertices_;
    typedef typename Value<TWaveletTreeVertices_>::Type TWaveletTreeVertex_;
    typedef typename Value<TWaveletTreeVertex_, 2>::Type TPos;

    typedef Pair<TChar, TPos> Type;
};

template <typename TChar, typename TSpec>
struct Value<RightArrayBinaryTree<TChar, TSpec> const>
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type TWaveletTreeVertices_;
    typedef typename Value<TWaveletTreeVertices_>::Type TWaveletTreeVertex_;
    typedef typename Value<TWaveletTreeVertex_, 2>::Type TPos;

    typedef Pair<TChar, TPos> Type;
};

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================

/**
.Class.RightArrayBinaryTree:
..cat:WaveletTree
..summary:A special format to encode the structure of a wavelet tree. The structure is very space efficient because only one position is stored which encodes where the left and right subtree of a given node exist.
..signature:RightArrayBinaryTree<TValue, TSpec>
..param.TSpec:The value type, that is the type of the stored characters.
..param.TSpec:The wavelet tree structure specialisation.
...default:void.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec>
class RightArrayBinaryTree
{
public:
    typename Fibre<RightArrayBinaryTree, FibreTreeStructureEncoding>::Type treeVertices;
    TChar minCharValue;
 
    RightArrayBinaryTree() :
        treeVertices(),
        minCharValue()
    {}

    template <typename TText>
    explicit RightArrayBinaryTree(TText const & text) :
        treeVertices(),
        minCharValue()
    {
        createRightArrayBinaryTree(*this,
                                    text);
    }

    inline bool operator==(const RightArrayBinaryTree & b) const
    {
        return treeVertices == b.treeVertices;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

/**
.Function.clear
..param.object:
...type:Class.RightArrayBinaryTree
*/
template <typename TChar, typename TSpec>
inline void clear(RightArrayBinaryTree<TChar, TSpec> & treeStructure)
{
    clear(treeStructure.treeVertices);
}

// ==========================================================================
/**
.Function.createRightArrayBinaryTree
..summary:Computes the wavelet tree structure of a text.
..signature:createRightArrayBinaryTree(waveletTreeStructure, text)
..param.waveletTreeStructure:A wavelet tree structure.
...type:Class.RightArrayBinaryTree
..param.text:A text.
...type:Class.String
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RightArrayBinaryTree<Dna5> waveletTreeStructure;
computeRightArrayBinaryTree(genome);
*/
// This function computes the wavelet tree structure.
template <typename TChar, typename TSpec, typename TBorderString, typename TPrefixSumTable, typename TIterSpec>
inline void createRightArrayBinaryTree(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it,
                                        TBorderString & borderString,
                                        TPrefixSumTable & pst)
{
    do
    {
        // TODO (singer): Comment this
        if (back(borderString).i2 - back(borderString).i1 + 1 < 3 ||
            getPrefixSum(pst, back(borderString).i1) == getPrefixSum(pst, back(borderString).i2 + 1))
        {
            setCharacter(it, getCharacter(pst, back(borderString).i1 + 1));
            SEQAN_ASSERT_MSG(isLeaf(it), "You just deleted a subtree.");
        }
        else
            _setChildVertices(it, borderString, pst);

        if (!_goDownConstruction(it) && !_setAndGoRight(it, borderString, pst))
            while (_goUpStructureConstruction(it, borderString) && !_setAndGoRight(it, borderString, pst))
                ;
    }
    while (!isRoot(it));
}

// ==========================================================================
// This function computes the wavelet tree structure.
template <typename TChar, typename TSpec, typename TIterSpec, typename TPrefixSumTable>
inline void createRightArrayBinaryTree(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it,
                                        TPrefixSumTable & pst)
{
    typedef RightArrayBinaryTree<TChar, TSpec> TRightArrayBinaryTree;
    TRightArrayBinaryTree & waveletTreeStructure = container(it);

    unsigned alpSize = getAlphabetSize(pst);
    String<Pair<unsigned> > borderString;
    appendValue(borderString, Pair<unsigned>(0, alpSize - 1));
    _resize(waveletTreeStructure, 1, Exact());
    createRightArrayBinaryTree(it, borderString, pst);
}

// This function computes the wavelet tree structure contained in the lfTable.
template < typename TValue, typename TSpec, typename TPrefixSumTable>
inline void createRightArrayBinaryTree(LfTable<SentinelRankDictionary<RankDictionary<WaveletTree<TValue> >, TSpec>, TPrefixSumTable> & lfTable)
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type TRightArrayBinaryTree;
    TRightArrayBinaryTree & rightArrayBinaryTree = lfTable.occTable.rankDictionary.waveletTreeStructure;

    typename Iterator<TRightArrayBinaryTree, TopDown<ParentLinks<void> > >::Type it(rightArrayBinaryTree, 0u);
    createRightArrayBinaryTree(it, lfTable.prefixSumTable);
}

template <typename TChar, typename TSpec, typename TText>
inline void createRightArrayBinaryTree(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TText const & text)
{
    PrefixSumTable<TChar, void> pst(text);

    typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<ParentLinks<> > >::Type it(waveletTreeStructure, 0u);
    createRightArrayBinaryTree(it, pst);
}

// ==========================================================================
/**
.Function.empty
..param.object:
...type:Class.RightArrayBinaryTree
*/
template <typename TChar, typename TSpec>
inline bool empty(RightArrayBinaryTree<TChar, TSpec> const & treeStructure)
{
    return empty(getFibre(treeStructure, FibreTreeStructureEncoding()));
}

// ==========================================================================
/**
.Function.RightArrayBinaryTree#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.RightArrayBinaryTree
..cat:Index
..param.container:The container holding the fibre.
...type:Class.RightArrayBinaryTree
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.RightArrayBinaryTree Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec>
inline typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type &
getFibre(RightArrayBinaryTree<TChar, TSpec> & treeStructure, FibreTreeStructureEncoding)
{
    return treeStructure.treeVertices;
}

template <typename TChar, typename TSpec>
inline typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type const &
getFibre(RightArrayBinaryTree<TChar, TSpec> const & treeStructure, FibreTreeStructureEncoding)
{
    return treeStructure.treeVertices;
}

// ==========================================================================
// This function returns the number of different entries in the wavelet tree structure.
template <typename TChar, typename TSpec>
inline unsigned _length(RightArrayBinaryTree<TChar, TSpec> const & tree)
{
    return length(tree.treeVertices);
}

// ==========================================================================
// This function resizes the string holding the nodes of the wavelet tree structure.
template <typename TChar, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RightArrayBinaryTree<TChar, TSpec> >::Type
_resize(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TSize size, Tag<TExpand> tag)
{
    return resize(treeStructure.treeVertices, size, tag);
}

// This function resizes the string holding the nodes of the wavelet tree structure.
template <typename TChar, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RightArrayBinaryTree<TChar, TSpec> >::Type
_resize(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TSize size,
                   typename Value<typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type>::Type value,
                   Tag<TExpand> tag)
{
    return resize(treeStructure.treeVertices, size, value, tag);
}

// ==========================================================================
/**
.Function.RightArrayBinaryTree#open
..class:Class.RightArrayBinaryTree
..summary:This functions loads a @Class.RightArrayBinaryTree@ from disk.
..signature:open(rightArrayBinaryTree, fileName [, openMode])
..param.rightArrayBinaryTree:The rightArrayBinaryTree.
...type:Class.RightArrayBinaryTree
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec>
inline bool open(
    RightArrayBinaryTree<TChar, TSpec> & treeStructure,
    const char * fileName,
    int openMode)
{
    String<char> name;

    String<TChar> minString;
    name = fileName;    append(name, ".rtv");   if(!open(getFibre(treeStructure, FibreTreeStructureEncoding()), toCString(name), openMode)) return false;
    name = fileName;    append(name, ".rtm");   if(!open(minString, toCString(name), openMode)) return false;
    treeStructure.minCharValue = minString[0];
    return true;
}

template <typename TChar, typename TSpec>
inline bool open(
    RightArrayBinaryTree<TChar, TSpec> & treeStructure,
    const char * fileName)
{
    return open(treeStructure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}

// ==========================================================================
/**
.Function.RightArrayBinaryTree#save
..class:Class.RightArrayBinaryTree
..summary:This functions saves a @Class.RightArrayBinaryTree@ to disk.
..signature:save(rightArrayBinaryTree, fileName [, openMode])
..param.rightArrayBinaryTree:The rightArrayBinaryTree.
...type:Class.RightArrayBinaryTree
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec>
inline bool save(
    RightArrayBinaryTree<TChar, TSpec> const & treeStructure,
    const char * fileName,
    int openMode)
{
    String<char> name;

    String<TChar> minString;
    appendValue(minString, treeStructure.minCharValue);
    name = fileName;    append(name, ".rtv");   if (!save(getFibre(treeStructure, FibreTreeStructureEncoding()), toCString(name), openMode)) return false;
    name = fileName;    append(name, ".rtm");   if (!save(minString, toCString(name), openMode)) return false;
    return true;
}

template <typename TChar, typename TSpec>
inline bool save(
    RightArrayBinaryTree<TChar, TSpec> const & treeStructure,
    const char * fileName)
{
    return save(treeStructure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}

}

#endif // INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H
