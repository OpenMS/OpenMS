// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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

#ifndef INDEX_FM_RANKDICTIONARY_WT
#define INDEX_FM_RANKDICTIONARY_WT

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

struct FibreBitStrings_;
struct FibreTreeStructure_;
struct FibreDollarPosition_;
struct FibreOccTable_;
struct FibreRankDictionary_;

template <typename TValue>
struct WaveletTree;

// ==========================================================================
// Tags
// ==========================================================================

typedef Tag<FibreTreeStructure_> const FibreTreeStructure;
typedef Tag<FibreBitStrings_> const FibreBitStrings;
typedef Tag<FibreDollarPosition_> const FibreDollarPosition;
typedef Tag<FibreOccTable_> const FibreOccTable;
typedef Tag<FibreRankDictionary_> const FibreRankDictionary;

// ==========================================================================
// Metafunctions
// ==========================================================================
/**
.Tag.WaveletTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Spec.WaveletTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a @Spec.WaveletTree@.
..cat:Spec.WaveletTree

..tag.FibreBitStrings:The string set containing a bit string for each node.

..tag.FibreTreeStructure:The wavelet tree structure of the wavelet tree.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

///.Metafunction.Fibre.param.TContainer.type:Class.RankDictionary
///.Metafunction.Fibre.param.TSpec.type:Tag.WaveletTree Fibres
/*
.Metafunction.Fibre:
..summary:Type of a specific FMIndex member (fibre).
..signature:Fibre<RankDictionary<WaveletTree<TValue> >, TFibreSpec>::Type
..class:Spec.FMIndex
..cat:Index
..param.TValue:The character type of the @Class.String@ the wavelet tree represents.
..param.TFibreSpec:Tag to specify the fibre.
...type:Tag.WaveletTree Fibres
..returns:Fibre type.
..remarks:Some containers, such as @Spec.FMIndex@, can be seen as a bundle consisting of various fibres. Because not 
every table is a fibre we did not call them tables, however, in many cases one can think of fibres as tables. The 
fibre interface was designed to unify the access to the members of the different fibres.
To get a reference or the type of a specific fibre use @Function.getFibre@ or @Metafunction.Fibre@.		
..include:seqan/index.h
*/

template <typename TValue>
struct Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>
{
    typedef StringSet<RankSupportBitString<void> > Type;
};

template <typename TValue>
struct Fibre<RankDictionary<WaveletTree<TValue> > const, FibreBitStrings>
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>::Type const Type;
};

template <typename TValue>
struct Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>
{
    typedef typename MakeUnsigned<TValue>::Type TUChar_;
    typedef RightArrayBinaryTree<TUChar_, void>  Type;
};

template <typename TValue>
struct Fibre<RankDictionary<WaveletTree<TValue> > const, FibreTreeStructure>
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------
 
template <typename TValue>
struct Size<RankDictionary<WaveletTree<TValue> > >
{
    typedef typename Size<String<TValue> >::Type Type;
};

template <typename TValue>
struct Size<RankDictionary<WaveletTree<TValue> > const> :
    public Size<RankDictionary<WaveletTree<TValue> > > {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue>
struct Value<RankDictionary<WaveletTree<TValue> > >
{
    typedef TValue Type;
};

template <typename TValue>
struct Value<RankDictionary<WaveletTree<TValue> > const> :
    public Value<RankDictionary<WaveletTree<TValue> > > {};


// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class RankDictionary
// ----------------------------------------------------------------------------

/**
.Class.RankDictionary:
..cat:Index
..summary:A rank dictionary is a data structure to store the rank of an element of a sequence at every position of the 
sequence.
..signature:RankDictionary<TSpec>
..param.TSpec:The rank dictionary specialisation.
...type:Spec.WaveletTree
...type:Spec.SequenceBitMask
...default:@Spec.WaveletTree@
..include:seqan/index.h
*/
template<typename TSpec> 
class RankDictionary;

// ----------------------------------------------------------------------------
// Spec WaveletTree
// ----------------------------------------------------------------------------

/**
.Spec.WaveletTree:
..cat:Index
..summary:A wavelet tree is a tree like binary encoding of a text.
..signature:WaveletTree<TValue>
..param.TValue:The value type of the wavelet tree.
..include:seqan/index.h
..remarks:The nodes of a wavelet tree consist of a bit string as well as a character c. In each level of the tree, 
characters smaller than c are represented as a 0 while character greater or equal to c are represented with a 1.
The characters represented by a 0 form the string to be represented by the left subtree while characters represented
by a 1 form the string of the right subtree. Therefore, only the bit string of the root node represents all characters while all other nodes represent subsets.
*/

template <typename TValue>
class RankDictionary<WaveletTree<TValue> >
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>::Type    TBitStrings;
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type TWaveletTreeStructure;

public:
    TBitStrings bitStrings;
    TWaveletTreeStructure waveletTreeStructure;

    RankDictionary() {}

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    }

    template <typename TText, typename TFreqTable>
    RankDictionary(TText const & text, TFreqTable const & freqTable)
    {
        createRankDictionary(*this, text, freqTable);
    }

    template <typename TText, typename TFreqTable, typename TPrefixSumTable>
    RankDictionary(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable)
    {
        createRankDictionary(*this, text, freqTable, prefixSumTable);
    }

    RankDictionary & operator=(RankDictionary const & other)
    {
        bitStrings = other.bitStrings;
        waveletTreeStructure = other.waveletTreeStructure;
        return *this;
    }

    bool operator==(RankDictionary const & b) const
    {
        typedef typename Size<typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>::Type>::Type TSize;
        
        if (length(bitStrings) != length(b.bitStrings))
            return false;

        for (TSize i = 0; i < length(bitStrings); ++i)
            if (!(bitStrings[i] == b.bitStrings[i]))
                return false;

        return waveletTreeStructure == b.waveletTreeStructure;
    }
};

// ==========================================================================
//Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#clear
..class:Class.RankDictionary
..summary:Clears the rank dictionary.
..signature:clear(dictionary)
..param.dictionary:The rank dictionary to be cleared.
...type:Class.RankDictionary
..include:seqan/index.h
*/

template <typename TValue>
inline void clear(RankDictionary<WaveletTree<TValue> > & dictionary)
{
    clear(getFibre(dictionary, FibreBitStrings()));
    clear(getFibre(dictionary, FibreTreeStructure()));
}

// ----------------------------------------------------------------------------
// Function empty
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#empty
..class:Class.RankDictionary
..summary:Returns whether or not the rank dictionary is empty.
..signature:empty(dictionary)
..param.dictionary:The rank dictionary to be checked.
...type:Class.RankDictionary
..returns:$true$ if the dictionary is empty, $false$ otherwise.
..include:seqan/index.h
*/

template <typename TValue>
inline bool empty(RankDictionary<WaveletTree<TValue> > const & dictionary)
{
    return empty(getFibre(dictionary, FibreTreeStructure()));
}

// ----------------------------------------------------------------------------
// Function getValue
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#getValue
..summary:Returns the character of a specified position.
..signature:getCharacter(dictionary, pos)
..class:Class.RankDictionary
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.pos:The position
..include:seqan/index.h
*/

template <typename TValue, typename TPos>
inline TValue
getValue(RankDictionary<WaveletTree<TValue> > & dictionary,
                 TPos pos)
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type const    TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type                 TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type                                       TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type                                     TChar;

    unsigned treePos = 0;
    typename Iterator<TWaveletTreeStructure, TopDown<> >::Type iter(dictionary.waveletTreeStructure, treePos);
    
    // initialize the return value with the smallest possible value
    TChar character = dictionary.waveletTreeStructure.minCharValue;

    // while the current node is not a leaf, go right if the bit at the current position is 1
    // go left otherwise
    // note that when going right the return value changes
    while (true)
    {
        TPos rank1 = getRank(dictionary.bitStrings[treePos], pos);
        if (isBitSet(dictionary.bitStrings[treePos], pos))
        {
            character = getCharacter(iter); 
            pos = rank1 - 1;  // -1 because strings start at 0
            if (!goRightChild(iter))
                break;
        }
        else
        {
            pos -= rank1;
            if (!goLeftChild(iter))
                break;
        }
        treePos = getPosition(iter);
    }

    return character;
}

// TODO(singer): We would like to have only one variant BUT there is a getValue() with const.
template <typename TValue, typename TPos>
inline TValue
getValue(RankDictionary<WaveletTree<TValue> > const & tree, TPos pos)
{
    return getValue(const_cast<RankDictionary<WaveletTree<TValue> > &>(tree), pos);
}

// ----------------------------------------------------------------------------
// Function getFibre
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#getFibre:
..summary:Returns a specific fibre of a dictionary.
..signature:getFibre(dictionary, fibreTag)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary holding the fibre.
...type:Class.RankDictionary
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.WaveletTree Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
template <typename TValue>
inline typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>::Type &
getFibre(RankDictionary<WaveletTree<TValue> >& dictionary, const FibreBitStrings)
{
    return dictionary.bitStrings;
}

template <typename TValue>
inline typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>::Type const &
getFibre(RankDictionary<WaveletTree<TValue> > const & dictionary, const FibreBitStrings)
{
    return dictionary.bitStrings;
}

template <typename TValue>
inline typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type &
getFibre(RankDictionary<WaveletTree<TValue> >& dictionary, FibreTreeStructure)
{
    return dictionary.waveletTreeStructure;
}

template <typename TValue>
inline typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type const &
getFibre(RankDictionary<WaveletTree<TValue> > const & dictionary, const FibreTreeStructure)
{
    return dictionary.waveletTreeStructure;
}

// ----------------------------------------------------------------------------
// Function countOccurrences
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#countOccurrences:
..summary:Returns the rank (number of occurrences) of a specified character up to a specified position. 
..signature:countOccurrences(dictionary, character, pos)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The character of interest.
..param.pos:The position (which is also included in the rank computation).
..returns:The rank (number of occurrences) of a specified character up to a specified position. 
...type:nolink:$unsigned$
..include:seqan/index.h
*/

template < typename TValue, typename TCharIn, typename TPos>
inline unsigned countOccurrences(RankDictionary<WaveletTree<TValue> > const & tree, TCharIn const character,
                                   TPos const pos)
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;

    TPos sum = pos;
    TPos treePos = 0;

    // determine the leaf containing the character
    // count the number of 1 or 0 up to the computed position
    typename Iterator<TWaveletTreeStructure const, TopDown<> >::Type it(tree.waveletTreeStructure, treePos);
    TChar charInTree = tree.waveletTreeStructure.minCharValue;
    while (true)
    {
        TPos addValue = getRank(tree.bitStrings[treePos], sum);
        if (ordGreater(getCharacter(it), character))
        {
            if (addValue > sum) return 0;

            sum -= addValue;
            if (!goLeftChild(it))
                break;
        }
        else
        {
            if (addValue == 0) return 0;

            charInTree = getCharacter(it);
            sum = addValue - 1;
            if (!goRightChild(it))
                break;
        }
        treePos = getPosition(it);
    }

    if (ordEqual(charInTree, character))
        return sum + 1;

    return 0;
}

// ----------------------------------------------------------------------------
// Function _fillWaveletTree
// ----------------------------------------------------------------------------

// This function is used to fill the bit strings of the wavelet tree.
template <typename TValue, typename TText>
inline void _fillWaveletTree(RankDictionary<WaveletTree<TValue> > & tree, TText const & text)
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>::Type     TFibreRankSupportBitStrings;
    typedef typename Value<TFibreRankSupportBitStrings>::Type                               TFibreRankSupportBitString;
    typedef typename Fibre<TFibreRankSupportBitString, FibreBits>::Type                     TFibreBitString;
    typedef typename Size<TFibreBitString>::Type                                            TSize;
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreTreeStructure>::Type  TWaveletTreeStructure;

    resize(tree.bitStrings, _length(tree.waveletTreeStructure), Exact());

    for (TSize i = 0; i < length(text); ++i)
    {
        typename Iterator<TWaveletTreeStructure, TopDown<> >::Type it(tree.waveletTreeStructure, 0);
        bool bit;

        while (true)
        {
            // decide whether the character is smaller then the pivot element of the current node 
            if (ordGreater(getCharacter(it), getValue(text, i)))
            {
                bit = 0;
                appendValue(getFibre(tree, FibreBitStrings())[getPosition(it)], bit);
                if (!goLeftChild(it))
                    break;
            }
            else
            {
                bit = 1;
                appendValue(getFibre(tree, FibreBitStrings())[getPosition(it)], bit);
                if (!goRightChild(it))
                    break;
            }
        }
    }

    TFibreRankSupportBitStrings & bitStrings = getFibre(tree, FibreBitStrings());
    for (TSize i = 0; i < length(bitStrings); ++i)
        _updateRanks(bitStrings[i]);
}

// ----------------------------------------------------------------------------
// Function createRankDictionary
// ----------------------------------------------------------------------------
/**
.Function.RankDictionary#createRankDictionary
..class:Class.RankDictionary
..summary:This functions creates the dictionary.
..signature:createRankDictionary(dictionary, text)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.text:A text to be transfered into a wavelet tree.
...type:Class.String
..include:seqan/index.h
*/

// TODO(singer): change to createRankDictionary
template <typename TValue, typename TText> 
inline void createRankDictionary(RankDictionary<WaveletTree<TValue> > & dictionary, TText const & text)
{
    createRightArrayBinaryTree(getFibre(dictionary, FibreTreeStructure()), text);
    _fillWaveletTree(dictionary, text);
}

template <typename TValue, typename TSpec, typename TPrefixSumTable, typename TText> 
inline void createRankDictionary(LfTable<SentinelRankDictionary<RankDictionary<WaveletTree<TValue> >, TSpec >, TPrefixSumTable> & lfTable,
                              TText const & text)
{
    createRightArrayBinaryTree(lfTable);
    _fillWaveletTree(getFibre(getFibre(lfTable, FibreOccTable()), FibreRankDictionary()), text);
}

// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#open
..class:Class.RankDictionary
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
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

template <typename TValue>
inline bool open(RankDictionary<WaveletTree<TValue> > & dictionary, const char * fileName, int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".wtc"); if (!open(getFibre(dictionary, FibreBitStrings()), toCString(name), openMode)) return false;
    name = fileName;    append(name, ".wts"); if (!open(getFibre(dictionary, FibreTreeStructure()), toCString(name), openMode)) return false;
    return true;
}

    template <typename TValue>
inline bool open(RankDictionary<WaveletTree<TValue> > & tree, const char * fileName)
{
    return open(tree, fileName, DefaultOpenMode<RankDictionary<WaveletTree<TValue> > >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#save
..class:Class.RankDictionary
..summary:This functions saves a dictionary to disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
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

template <typename TValue>
inline bool save(RankDictionary<WaveletTree<TValue> > const & dictionary, const char * fileName, int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".wtc");   
    if (!save(getFibre(dictionary, FibreBitStrings()), toCString(name), openMode)) return false;
    
    name = fileName;    append(name, ".wts");   
    if (!save(getFibre(dictionary, FibreTreeStructure()), toCString(name), openMode)) return false;
    
    return true;
}

template <typename TValue>
inline bool save(RankDictionary<WaveletTree<TValue> > const & tree, const char * fileName)
{
    return save(tree, fileName, DefaultOpenMode<RankDictionary<WaveletTree<TValue> > >::VALUE);
}

}
#endif  // INDEX_FM_RANKDICTIONARY_WT
