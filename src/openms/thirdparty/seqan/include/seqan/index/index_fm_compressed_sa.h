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


#ifndef INDEX_FM_COMPRESSED_SA_H_
#define INDEX_FM_COMPRESSED_SA_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

/**
.Tag.CompressedSA Fibres
..summary:Tag to select a specific fibre of a @Class.CompressedSA@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a sparse string.
..cat:Index

..tag.FibreSparseString:The sparse string.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
struct FibreSparseString_;
typedef Tag<FibreSparseString_> const FibreSparseString;

template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA;

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>
{
    typedef TSparseString Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Fibre<CompressedSA<TSparseString, TLfTable, TSpec> const, FibreSparseString>
{
    typedef TSparseString Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<CompressedSA<TSparseString, TLfTable, TSpec> >
{
    // TODO(singer): We actually need a proxy here.
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<const CompressedSA<TSparseString, TLfTable, TSpec> >
{
    // TODO(singer): We actually need a proxy here.
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type /*const*/ Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Value<CompressedSA<TSparseString, TLfTable, TSpec> >
{
    typedef typename Value<TSparseString>::Type Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Value<CompressedSA<TSparseString, TLfTable, TSpec> const>
{
    typedef typename Value<TSparseString>::Type const Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// forwards
template <typename TPos, typename TOffSet>
TPos _addGapDistance(TPos const & value, TOffSet const & offSet);

template <typename TSeqId, typename TSpec, typename TPos, typename TOffSet>
Pair<TSeqId, TPos> _addGapDistance(Pair<TSeqId, TPos, TSpec> const & value, TOffSet const & offSet);

// ----------------------------------------------------------------------------
// Class CompressedSA
// ----------------------------------------------------------------------------

/**
.Class.CompressedSA:
..cat:Index
..summary:A suffix array storing only a few suffix array entries and computing the remaining on demand.
..signature:CompressedSA<TSparseString, TLfTable, TSpec>
..param.TSparseString:The string containing specific suffix array entries.
...type:Class.SparseString
..param.TLfTable:The lfTable containg an occurrence table and a prefix sum table.
...type:Class.LfTable
..param.TSpec:Possibility to specialise a compressed suffix array.
...default:void.
..remarks:The compressed suffix array can only be used with the FM index.
..include:seqan/index.h
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA
{
public:
    TSparseString   sparseString;
    TLfTable *      lfTable;

    CompressedSA(){};

    CompressedSA(TLfTable & lfTable) :
        lfTable(&lfTable)
    {}

    inline CompressedSA & operator=(CompressedSA const & other)
    {
        sparseString = other.sparseString;
        lfTable = other.lfTable;
        return *this;
    }

    typedef typename Value<typename Fibre<TSparseString, FibreValueString>::Type>::Type TCompressedSaValue;
    typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;

    template <typename TPos>
    inline TCompressedSaValue const operator[](TPos pos)
    {
        TIndicatorString const & indicatorString = getFibre(sparseString, FibreIndicatorString());
        TPos counter = 0;

        while (!isBitSet(indicatorString, pos))
        {
            pos = lfMapping(*lfTable, pos);
            ++counter;
        }
        return _addGapDistance(getValue(sparseString.valueString, getRank(indicatorString, pos) - 1), counter);
    }

    template <typename TPos>
    inline TCompressedSaValue operator[](TPos pos) const
    {
        TIndicatorString const & indicatorString = getFibre(sparseString, FibreIndicatorString());
        TPos counter = 0;
        while (!isBitSet(indicatorString, pos))
        {
            pos = lfMapping(*lfTable, pos);
            ++counter;
        }
        return _addGapDistance(getValue(sparseString.valueString, getRank(indicatorString, pos) - 1), counter);
    }

    inline bool operator==(const CompressedSA & other) const
    {
        return sparseString == other.sparseString &&
               *lfTable == *(other.lfTable);
    }

};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function _addGapDistance
// ----------------------------------------------------------------------------

template <typename TPos, typename TOffSet>
TPos _addGapDistance(TPos const & value, TOffSet const & offSet)
{
    return value + offSet;
}

template <typename TSeqId, typename TSpec, typename TPos, typename TOffSet>
Pair<TSeqId, TPos> _addGapDistance(Pair<TSeqId, TPos, TSpec> const & value, TOffSet const & offSet)
{
    return Pair<TSeqId, TPos>(value.i1, value.i2 + offSet);
}

// ==========================================================================
/*
.Function.assignValue
..param.object
...type:Class.CompressedSA
..remarks:In the case of the compressed suffix array the new value only influences the sparse string and not the indicator string.
..include:seqan/index.h
*/
// template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos, typename TValue>
// void assignValue(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos, TValue value)
// {
//     assignValue(getFibre(compressedSA, FibreSparseString()), pos, value);
// }

// ==========================================================================
/**
.Function.clear.param.object.type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline void clear(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    clear(getFibre(compressedSA, FibreSparseString()));
}

// ==========================================================================
/**
.Function.empty.param.object.type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool empty(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return empty(getFibre(compressedSA, FibreSparseString()));
}

// ==========================================================================
// TODO (singer): How to name this function?
/*
.Function.entryStored
..summary:Determines whether the requested position contains a value different from the default value.
..signature:bool entryStored(object, pos)
..param.object
...type:Class.CompressedSA
..param.pos:The position of an item in object.
...Remarks:pos should be convertible to Position<T>::Type for container-type T.
..include:seqan/index.h
*/
template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline bool entryStored(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, TPos const & pos)
{
    return entryStored(getFibre(compressedSA, FibreSparseString()), pos);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline bool entryStored(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos const & pos)
{
    return entryStored(getFibre(compressedSA, FibreSparseString()), pos);
}

// ==========================================================================
// This function creates a compressed suffix array using a normal one.
/**
.Function.CompressedSA#createCompressedSa
..summary:This functions creates a compressed suffix array with a specified compression factor.
..signature:void createCompressedSa(compressedSA, completeSA, compressionFactor [,offset])
..param.compressedSA:The compressed suffix array
...type:Class.CompressedSA
..param.completeSA:A complete suffix array containing all values
..param.compressionFactor:The compression factor.
...type:Concept.UnsignedIntegerConcept
...remarks:A compression factor of x means that the compressed suffix array specifically stores a value for every x values in the complete suffix array.
..param:offset:Number of elements at the beginning which should contain the default value.
..include:seqan/index.h
*/
template <typename TSparseString, typename TLfTable, typename TSpec, typename TSA, typename TCompression, typename TSize>
void createCompressedSa(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, 
        TSA const & sa,
        TCompression const compressionFactor, 
        TSize offset)
{
    typedef CompressedSA<TSparseString, TLfTable, TSpec>            TCompressedSA;
    typedef typename Size<TSA>::Type                                TSASize;
    typedef typename Fibre<TCompressedSA, FibreSparseString>::Type  TSparseSA;
    typedef typename Fibre<TSparseSA, FibreIndicatorString>::Type   TIndicatorString;
    typedef typename Iterator<TSA const, Standard>::Type            TSAIter;

    TSparseSA & sparseString = getFibre(compressedSA, FibreSparseString());
    TIndicatorString & indicatorString = getFibre(sparseString, FibreIndicatorString());

    TSASize saLen = length(sa);
    resize(compressedSA, saLen + offset, Exact());
    
    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());

    for (TSASize pos = offset; saIt != saItEnd; ++saIt, ++pos)
    {
        if (getSeqOffset(getValue(saIt)) % compressionFactor == 0)
            setBit(indicatorString, pos);
        else
            clearBit(indicatorString, pos);
    }
    _updateRanks(indicatorString);

    resize(sparseString.valueString, getRank(indicatorString, length(indicatorString) - 1), Exact());

    saIt = begin(sa, Standard());
    for (TSASize pos = offset, counter = 0; saIt != saItEnd; ++saIt, ++pos)
    {
        if (isBitSet(indicatorString, pos))
        {
            assignValue(compressedSA.sparseString.valueString, counter, getValue(saIt));
            ++counter;
        }
    }
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSA, typename TCompression>
void createCompressedSa(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TSA const & completeSA, TCompression const compressionFactor)
{
    createCompressedSa(compressedSA, completeSA, compressionFactor, 0);
}

// ==========================================================================
/**
.Function.CompressedSA#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.CompressedSA
..cat:Index
..param.container:The container holding the fibre.
...type:Class.CompressedSA
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.CompressedSA Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type const &
getFibre(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, FibreSparseString)
{
    return compressedSA.sparseString;
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type &
getFibre(CompressedSA<TSparseString, TLfTable, TSpec>&compressedSA, FibreSparseString)
{
    return compressedSA.sparseString;
}

// ==========================================================================
// This functions computes the position in the suffix array of text[sa[pos] - 1]
// iff the current position is not present in the compressed suffix array.
template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline bool _getNextPos(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, TPos & pos)
{
    typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;
    TIndicatorString const & indicatorString = compressedSA.sparseString.indicatorString;

    if (isBitSet(indicatorString, pos))
    {
        return true;
    }
    pos = lfMapping(*compressedSA.lfTable, pos);
    return false;
}

// ==========================================================================
/**
.Function.length.param.object.type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Size<typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type>::Type
length(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA)
{
    return length(getFibre(compressedSA, FibreSparseString()));
}

// ==========================================================================
/**
.Function.resize.param.object.type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize, typename TExpand>
inline typename Size<typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type>::Type
resize(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
                   TSize size,
                   Tag<TExpand> tag)
{
    return resize(getFibre(compressedSA, FibreSparseString()), size, tag);
}

// template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
// inline void resize(CompressedSA<String<TSparseString>, TLfTable, TSpec> & compressedSA, TSize size)
// {
//     resize(getFibre(compressedSA, FibreSparseString()), size);
// }

// ==========================================================================
/**
.Function.open
..param.object:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool open(
    CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;
    if (!open(getFibre(compressedSA, FibreSparseString()), toCString(name), openMode))
    {
        return false;
    }
    return true;

}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool open(
    CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
    const char * fileName)
{
    return open(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
}

// ==========================================================================
/**
.Function.CompressedSA#save
..class:Class.CompressedSA
..summary:This functions saves a compressed suffix array to disk.
..signature:open(compressedSA, fileName [, openMode])
..param.compressedSA:The string to be saved.
...type:Class.CompressedSA
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
template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool save(
    CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;
    if (!save(getFibre(compressedSA, FibreSparseString()), toCString(name), openMode))
    {
        return false;
    }
    return true;
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool save(
    CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
    const char * fileName)
{
    return save(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
}

// ==========================================================================
/**
.Function.setLfTable
..summary:Set the LfTable of the compressed suffix array.
..signature:setLfTable(CompressedSA<TSparseString, TLfTable, TSpec> compressedSa, TLfTable & lfTable)
..param.CompressedSA<TSparseString, TLfTable, TSpec>:The compressed suffix array.
...type:Class.CompressedSA
..param.lfTable
...type:Class.LfTable
..include:seqan/index.h
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
void setLfTable(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TLfTable & lfTable)
{
    compressedSA.lfTable = &lfTable;
}

// ==========================================================================
/**
.Function.CompressedSA#value
..summary:Returns the value stored at a specified position in the compressed suffix-array.
..signature:value(compressedSA, pos)
..param.compressedSA:The compressed suffix array to access.
...type:Class.CompressedSA
..param.pos:Position at which to access the suffix array.
...type:Concept.UnsignedIntegerConcept
..remarks:Note that the compressed suffix array is read only. Therefore a const reference is return by
this function.
*/
template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<TSparseString>::Type
value(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
    return compressedSA[pos];
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<TSparseString>::Type const
value(const CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
    return compressedSA[pos];
}

}

#endif // INDEX_FM_COMPRESSED_SA_H_
