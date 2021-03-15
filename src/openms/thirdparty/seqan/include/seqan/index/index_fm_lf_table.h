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

#ifndef INDEX_FM_LF_TABLE_H_
#define INDEX_FM_LF_TABLE_H_

namespace seqan {

// ==========================================================================
// LfTable is an object storing all necessary information for the LF-mapping.
// To be more precise, the occurrence-table data structure as well as the
// prefix-sum table are stored.
// ==========================================================================

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TValue>
struct WaveletTree;

template<typename TValue> 
class RankDictionary;

template<typename TRankDirectorySpec, typename TSpec> 
class SentinelRankDictionary;

/**
.Tag.LF Table Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibreOccTable:The occurrence table of the lf table.
..tag.FMTablePrefixSumTable:The prefix sum table of the lf table.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/
struct FibreOccTable_;
struct FibrePrefixSumTable_;

typedef Tag<FibreOccTable_> const FibreOccTable;
typedef Tag<FibrePrefixSumTable_> const FMTablePrefixSumTable;

template <typename TOccTable, typename TPrefixSumTable>
struct LfTable;

// ==========================================================================
// Tags
// ==========================================================================

struct FibreRankDictionary_;
typedef Tag<FibreRankDictionary_> const FibreRankDictionary;

struct FibreSentinentalPosition_;
typedef Tag<FibreSentinentalPosition_> const FibreSentinentalPosition;

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable>, FibreOccTable>
{
    typedef TOccTable Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable> const, FibreOccTable>
{
    typedef TOccTable const Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>
{
    typedef TPrefixSumTable Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable> const, FMTablePrefixSumTable>
{
    typedef TPrefixSumTable const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TOccTable, typename TPrefixSumTable>
struct Reference<LfTable<TOccTable, TPrefixSumTable> >
{
    typedef TPrefixSumTable & Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class LfTable
// ----------------------------------------------------------------------------

/**
.Class.LfTable:
..cat:Index
..summary:LfTable is an object storing all necessary information for the LF-mapping.
..signature:LfTable<TOccTable, TPrefixSumTable>
..param.TOccTable:The occurrence table data structure.
...type:Spec.WaveletTree
..param.TPrefixSumTable:The specialisation tag.
...default:String
..include:seqan/Index.h
*/
template <typename TOccTable, typename TPrefixSumTable>
struct LfTable
{
    TOccTable occTable;
    TPrefixSumTable prefixSumTable;

    LfTable() :
        occTable(),
        prefixSumTable()
    {}

    LfTable(TOccTable const & occTable, TPrefixSumTable const & prefixSumTable) :
        occTable(occTable),
        prefixSumTable(prefixSumTable)
    {}

    inline LfTable & operator=(LfTable const & other)
    {
        occTable = other.occTable;
        prefixSumTable = other.prefixSumTable;
        return *this;
    }

    inline bool operator==(const LfTable & b) const
    {
        return occTable == b.occTable &&
               prefixSumTable == b.prefixSumTable;
    }

};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

/**
.Function.LfTable#clear
..class:Class.LfTable
..summary:Clears the LF table.
..signature:clear(lfTable)
..param.lfTable:The LF table to be cleared.
...type:Class.LfTable
..include:seqan/index.h
*/

template <typename TOccTable, typename TPrefixSumTable>
inline void clear(LfTable<TOccTable, TPrefixSumTable> & lfTable)
{
    clear(lfTable.occTable);
    clear(lfTable.prefixSumTable);
}

// ----------------------------------------------------------------------------
// Function empty
// ----------------------------------------------------------------------------

/**
.Function.LfTable#empty
..class:Class.LfTable
..summary:Clears the LF table.
..signature:empty(lfTable)
..param.lfTable:The LF table to be cleared.
...type:Class.LfTable
..returns:$true$ if the LF table is empty, $false$ otherwise.
...type:nolink:$bool$
..include:seqan/index.h
*/

template <typename TOccTable, typename TPrefixSumTable>
inline bool empty(LfTable<TOccTable, TPrefixSumTable> & lfTable)
{
    return empty(lfTable.occTable)
           && empty(lfTable.prefixSumTable);
}

// ----------------------------------------------------------------------------
// Function createLfTable
// ----------------------------------------------------------------------------

/**
.Function.createLfTable
..summary:Creates the LF table
..signature:createLfTable(lfTable, text)
..param.lfTable:The LF table to be constructed.
...type:Class.LfTable.
..param.text:The underlying text
...type:Class.String
..returns:$true$ on successes, $false$ otherwise.
..include:seqan/index.h
*/

template <typename TValue, typename TSpec, typename TPrefixSumTable, typename TText>
inline bool createLfTable(LfTable<SentinelRankDictionary<RankDictionary<WaveletTree<TValue> >, TSpec>, TPrefixSumTable> & lfTable, TText const text)
{
    typedef typename SAValue<TText>::Type   TSAValue;
    typedef TValue                          TAlphabet;
    typedef typename Fibre<SentinelRankDictionary<RankDictionary<WaveletTree<TValue> >, TSpec>, FibreSentinentalPosition>::Type TDollarPos;

    String<TSAValue> sa;
    resize(sa, length(text), Exact());
    createSuffixArray(sa, text, Skew7());

    createPrefixSumTable(lfTable.prefixSumTable, text);

    TAlphabet dollarSub;
    _determineDollarSubstitute(lfTable.prefixSumTable, dollarSub);

    String<TAlphabet> bwt;
    TDollarPos dollarPos = 0;
    resize(bwt, _computeBwtLength(text), Exact());
    _createBwTable(bwt, dollarPos, text, sa, dollarSub);

    clear(sa);

    createOccurrenceTable(lfTable, bwt, dollarSub, dollarPos);
    clear(bwt);

    _insertDollar(lfTable.prefixSumTable, countSequences(text));

    //here comes the dollar modification
    //addDollarNode(index.lfTable.occTable, dollarSub, dollarPos);
    return true;
}

// ----------------------------------------------------------------------------
// Function getFibre
// ----------------------------------------------------------------------------

/**
.Function.LfTable#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.LfTable
..cat:Index
..param.container:The container holding the fibre.
...type:Class.LfTable
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.LF Table Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>::Type &
getFibre(LfTable<TOccTable, TPrefixSumTable>&lfTable, FMTablePrefixSumTable)
{
    return lfTable.prefixSumTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>::Type const &
getFibre(LfTable<TOccTable, TPrefixSumTable> const & lfTable, FMTablePrefixSumTable)
{
    return lfTable.prefixSumTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type &
getFibre(LfTable<TOccTable, TPrefixSumTable>&lfTable, FibreOccTable)
{
    return lfTable.occTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type const &
getFibre(LfTable<TOccTable, TPrefixSumTable> const & lfTable, FibreOccTable )
{
    return lfTable.occTable;
}

// ----------------------------------------------------------------------------
// Function lfMapping
// ----------------------------------------------------------------------------

/**
.Function.LfTable#lfMapping:
..summary:Returns the position of an character at a specified position of L in F. L corresponds to the last column of 
the sorted cyclic rotations of the original text, while F correspond to the first column.
..cat:Index
..signature:lfMapping(lfTable, pos)
..param.lfTable:The @Class.LfTable@ holding the occurrence and prefix sum table.
...type:Class.LfTable
..param.pos:The position in L
..returns:Returns the position of the character L[c] in F. The returned position is of the same type as pos.
..include:seqan/index.h
*/

template <typename TLfTable, typename TPos>
inline TPos lfMapping(TLfTable & lfTable,
                      TPos pos)
{
    typedef typename Fibre<TLfTable, FibreOccTable>::Type TOccTable;
    typedef typename Value<TOccTable>::Type TChar;
    TChar c = getValue(lfTable.occTable, pos);
    return countOccurrences(getFibre(lfTable, FibreOccTable()), c, pos) + getPrefixSum(lfTable.prefixSumTable, getCharacterPosition(lfTable.prefixSumTable, c)) - 1;
}

// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------

/**
.Function.LfTable#open
..class:Class.LfTable
..summary:This functions loads a LF table from disk.
..signature:open(lfTable, fileName [, openMode])
..param.lfTable:The lfTable.
...type:Class.LfTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A nolink:$bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TOccTable, typename TPrefixSumTable>
inline bool open(
    LfTable<TOccTable, TPrefixSumTable> & lfTable,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName; 
    if (!open(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
    {
        return false;
    }
    name = fileName;    open(getFibre(lfTable, FMTablePrefixSumTable()), toCString(name), openMode);
    return true;

}

template <typename TOccTable, typename TPrefixSumTable>
inline bool open(
    LfTable<TOccTable, TPrefixSumTable> & lfTable,
    const char * fileName)
{
    return open(lfTable, fileName, DefaultOpenMode<LfTable<TOccTable, TPrefixSumTable> >::VALUE);
}

/**
.Function.LfTable#save
..class:Class.LfTable
..summary:This functions saves a LF table to disk.
..signature:save(lfTable, fileName [, openMode])
..param.lfTable:The dictionary.
...type:Class.LfTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A nolink:$bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TOccTable, typename TPrefixSumTable>
inline bool save(
    LfTable<TOccTable, TPrefixSumTable> const & lfTable,
    const char * fileName,
    int openMode)
{
    String<char> name;

    name = fileName;
    if (!save(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
    {
        return false;
    }
    name = fileName;    save(getFibre(lfTable, FMTablePrefixSumTable()), toCString(name), openMode);
    return true;
}

template <typename TOccTable, typename TPrefixSumTable>
inline bool save(
    LfTable<TOccTable, TPrefixSumTable> const & lfTable,
    const char * fileName)
{
    return save(lfTable, fileName, DefaultOpenMode<LfTable<TOccTable, TPrefixSumTable> >::VALUE);
}

}
#endif // INDEX_FM_LF_TABLE_H_
