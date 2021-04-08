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

#ifndef INDEX_FM_SPARSE_STRING_H_
#define INDEX_FM_SPARSE_STRING_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TString, typename TSpec>
struct SparseString;

struct FibreValueString_;
struct FibreIndicatorString_;

typedef Tag<FibreValueString_>       const FibreValueString;
typedef Tag<FibreIndicatorString_>    const FibreIndicatorString;

// ==========================================================================
// Metafunctions
// ==========================================================================

template <typename TFibreValueString, typename TSpec>
struct Value<SparseString<TFibreValueString, TSpec> >
{
    typedef typename Value<TFibreValueString>::Type Type;
};

template <typename TFibreValueString, typename TSpec>
struct Value<SparseString<TFibreValueString, TSpec> const>
{
    typedef typename Value<TFibreValueString>::Type Type;
};

// ==========================================================================
template <typename TFibreValueString, typename TSpec>
struct GetValue<SparseString<TFibreValueString, TSpec> > :
    Value<SparseString<TFibreValueString, TSpec> > {};

template <typename TFibreValueString, typename TSpec>
struct GetValue<SparseString<TFibreValueString, TSpec> const> :
    Value<SparseString<TFibreValueString, TSpec> const> {};

// ==========================================================================
template <typename TFibreValueString, typename TSpec>
struct Reference<SparseString<TFibreValueString, TSpec> >
{
    typedef typename Value<SparseString<TFibreValueString, TSpec> >::Type Type;
};

template <typename TFibreValueString, typename TSpec>
struct Reference<SparseString<TFibreValueString, TSpec> const>
{
    typedef typename Value<SparseString<TFibreValueString, TSpec> >::Type const Type;
};

// ==========================================================================
template <typename TSpec>
struct DefaultValue {};

template <typename TFibreValueString, typename TSpec>
struct DefaultValue<SparseString<TFibreValueString, TSpec> const>
{
    typedef typename GetValue<SparseString<TFibreValueString, TSpec> const>::Type Type;
    static const Type VALUE = -1;
};

template <typename TFibreValueString, typename TSpec>
struct DefaultValue<SparseString<TFibreValueString, TSpec> >
{
    typedef typename GetValue<SparseString<TFibreValueString, TSpec> const>::Type Type;
    static const Type VALUE = -1;
};

// ==========================================================================
template <typename TFibreValueString, typename TSpec>
struct Fibre<SparseString<TFibreValueString, TSpec>, FibreValueString>
{
    typedef TFibreValueString Type;
};

template <typename TFibreValueString, typename TSpec>
struct Fibre<SparseString<TFibreValueString, TSpec> const, FibreValueString>
{
    typedef TFibreValueString const Type;
};

template <typename TFibreValueString, typename TSpec>
struct Fibre<SparseString<TFibreValueString, TSpec>, FibreIndicatorString>
{
    typedef RankSupportBitString<void> Type;
};

template <typename TFibreValueString, typename TSpec>
struct Fibre<SparseString<TFibreValueString, TSpec> const, FibreIndicatorString>
{
    typedef RankSupportBitString<void> const Type;
};

// ==========================================================================
template <typename TFibreValueString, typename TSpec>
struct Iterator<SparseString<TFibreValueString, TSpec> const, Standard>
{
    typedef Iter<SparseString<TFibreValueString, TSpec> const, PositionIterator> Type;
};

template <typename TFibreValueString, typename TSpec>
struct Iterator<SparseString<TFibreValueString, TSpec>, Standard>
{
    typedef Iter<SparseString<TFibreValueString, TSpec>, PositionIterator> Type;
};

template <typename TFibreValueString, typename TSpec>
struct Iterator<SparseString<TFibreValueString, TSpec>, Rooted>:
    Iterator<SparseString<TFibreValueString, TSpec>, Standard>{};

template <typename TFibreValueString, typename TSpec>
struct Iterator<SparseString<TFibreValueString, TSpec> const, Rooted>:
    Iterator<SparseString<TFibreValueString, TSpec> const, Standard>{};



// ==========================================================================
// Classes
// ==========================================================================

/**
.Class.SparseString:
..cat:Index
..summary:A string storing only a fraction of the values of the original string..
..signature:SparseString<TValueString, TSpec>
..param.TValueString:The string containing the values.
..param.TSpec:The specialisation tag.
...default:void.
..include:seqan/String.h
*/
template <typename TValueString, typename TSpec = void>
struct SparseString
{
    typedef typename Fibre<SparseString, FibreValueString>::Type TFibreValueString_;
    typedef typename Fibre<SparseString, FibreIndicatorString>::Type TFibreIndicatorString;

    TFibreValueString_              	valueString;
    TFibreIndicatorString               indicatorString;

    inline bool operator==(const SparseString & b) const
    {
        return valueString == b.valueString &&
            indicatorString == b.indicatorString;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

/*
.Function._assignValueInValueString
..param.container:
...type:Class.CompressedSA
*/
template <typename TFibreValueString, typename TSpec, typename TPos, typename TValue>
inline void _assignValueInValueString(SparseString<TFibreValueString, TSpec> & string, TPos pos, TValue value)
{
    getFibre(string, FibreValueString())[pos] = value;
}

///.Function.clear.param.object.type:Class.SparseString
template <typename TFibreValueString, typename TSpec>
inline void clear(SparseString<TFibreValueString, TSpec> & string)
{
    clear(getFibre(string, FibreValueString()));
    clear(getFibre(string, FibreIndicatorString()));
}

///.Function.empty.param.object.type:Class.SparseString
template <typename TFibreValueString, typename TSpec>
inline bool empty(SparseString<TFibreValueString, TSpec> const & string)
{
    return empty(getFibre(string, FibreIndicatorString()));
}

template <typename TFibreValueString, typename TSpec, typename TPos>
inline bool _isContained(SparseString<TFibreValueString, TSpec> const & string, TPos const & pos)
{
    return isBitSet(getFibre(string, FibreIndicatorString()), pos);
}

///.Function.assignValue.param.container.type:Class.SparseString
template <typename TFibreValueString, typename TSpec, typename TPos, typename TValue>
inline void
assignValue(SparseString<TFibreValueString, TSpec> & string, TPos pos, TValue value)
{
    if (!_isContained(string, pos))
    {
        setBit(getFibre(string, FibreIndicatorString()), pos, 0);    
    }
    getFibre(string, FibreValueString())[getRank(getFibre(string, FibreIndicatorString()), pos) - 1] = value;
}

///.Function.getValue.param.container.type:Class.SparseString
template <typename TFibreValueString, typename TSpec, typename TPos>
inline typename GetValue<SparseString<TFibreValueString, TSpec> >::Type
getValue(SparseString<TFibreValueString, TSpec> & string, TPos pos)
{
    if (_isContained(string, pos))
        return getValue(getFibre(string, FibreValueString()), getRank(getFibre(string, FibreIndicatorString()), pos) - 1);
    else 
        return DefaultValue<SparseString<TFibreValueString, TSpec> >::VALUE;
}

template <typename TFibreValueString, typename TSpec, typename TPos>
inline typename GetValue<SparseString<TFibreValueString, TSpec> const>::Type
getValue(SparseString<TFibreValueString, TSpec> const & string, TPos pos)
{
    if (_isContained(string, pos))
        return getValue(getFibre(string, FibreValueString()), getRank(getFibre(string, FibreIndicatorString()), pos) - 1);
    else 
        return DefaultValue<SparseString<TFibreValueString, TSpec> const>::VALUE;
}


///.Function.value.param.container.type:Class.SparseString
template <typename TFibreValueString, typename TSpec, typename TPos>
inline typename Reference<SparseString<TFibreValueString, TSpec> >::Type 
value(SparseString<TFibreValueString, TSpec>&string, TPos pos)
{
    return getValue(string, pos);
}

template <typename TFibreValueString, typename TSpec, typename TPos>
inline typename Reference<SparseString<TFibreValueString, TSpec> >::Type
value(SparseString<TFibreValueString, TSpec> const & string, TPos pos)
{
    return getValue(string, pos);
}

///.Function.getFibre.param.container.type:Class.CompressedSA
template <typename TFibreValueString, typename TSpec>
inline typename Fibre<SparseString<TFibreValueString, TSpec>, FibreValueString>::Type const &
getFibre(SparseString<TFibreValueString, TSpec> const & sparseString, FibreValueString)
{
    return sparseString.valueString;
}

template <typename TFibreValueString, typename TSpec>
inline typename Fibre<SparseString<TFibreValueString, TSpec>, FibreValueString>::Type &
getFibre(SparseString<TFibreValueString, TSpec>&sparseString, FibreValueString)
{
    return sparseString.valueString;
}

template <typename TFibreValueString, typename TSpec>
inline typename Fibre<SparseString<TFibreValueString, TSpec>, FibreIndicatorString>::Type const &
getFibre(SparseString<TFibreValueString, TSpec> const & sparseString, FibreIndicatorString)
{
    return sparseString.indicatorString;
}

template <typename TFibreValueString, typename TSpec>
inline typename Fibre<SparseString<TFibreValueString, TSpec>, FibreIndicatorString>::Type &
getFibre(SparseString<TFibreValueString, TSpec>&sparseString, FibreIndicatorString)
{
    return sparseString.indicatorString;
}

///.Function.length.param.object.type:Class.SparseString
template <typename TFibreValueString, typename TSpec>
inline typename Size<typename Fibre<SparseString<TFibreValueString, TSpec>, FibreValueString>::Type>::Type
length(SparseString<TFibreValueString, TSpec> const & string)
{
    return length(getFibre(string, FibreIndicatorString()));
}

///.Function.resize.param.object.type:Class.SparseString
template <typename TFibreValueString, typename TSpec, typename TSize, typename TValue, typename TExpand>
inline typename Size<typename Fibre<SparseString<TFibreValueString, TSpec>, FibreValueString>::Type>::Type
resize(SparseString<TFibreValueString, TSpec> & string,
                   TSize const size,
                   TValue const value,
                   Tag<TExpand> const tag)
{
    if (value != DefaultValue<SparseString<TFibreValueString, TSpec> >::VALUE)
    {
        TSize _length = length(getFibre(string, FibreIndicatorString()));
        if (_length < size)
            resize(getFibre(string, FibreValueString()), length(getFibre(string, FibreValueString())) + (size - length), value, tag);
        else
            resize(getFibre(string, FibreValueString()), getRank(string, size), value, tag);
        resize(getFibre(string, FibreIndicatorString()), size, 1, tag);            
    }
    return resize(getFibre(string, FibreIndicatorString()), size, 0);
}

template <typename TFibreValueString, typename TSpec, typename TSize, typename TExpand>
inline typename Size<typename Fibre<SparseString<TFibreValueString, TSpec>, FibreValueString>::Type>::Type
resize(SparseString<TFibreValueString, TSpec> & string,
                   TSize const size,
                   Tag<TExpand> tag)
{
    return resize(getFibre(string, FibreIndicatorString()), size, 0, tag);
}


/**
.Function.open
..param.string:
...type:Class.SparseString
*/
template <typename TFibreValueString, typename TSpec>
inline bool open(
    SparseString<TFibreValueString, TSpec> & sparseString,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".val");
    if (!open(getFibre(sparseString, FibreValueString()), toCString(name), openMode)) // val = value string
    {
        return false;
    }
    name = fileName;    append(name, ".ind");   open(getFibre(sparseString, FibreIndicatorString()), toCString(name), openMode); // ind = indicator string
    return true;
}

template <typename TFibreValueString, typename TSpec>
inline bool open(
    SparseString<TFibreValueString, TSpec> & sparseString,
    const char * fileName)
{
    return open(sparseString, fileName, DefaultOpenMode<SparseString<TFibreValueString, TSpec> >::VALUE);
}

/**
.Function.SparseString#save
..class:Class.SparseString
..summary:This functions saves a sparse string to disk.
..signature:open(string, fileName [, openMode])
..param.string:The string to be saved.
...type:Class.SparseString
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
template <typename TFibreValueString, typename TSpec>
inline bool save(
    SparseString<TFibreValueString, TSpec> const & sparseString,
    const char * fileName)
{
    return save(sparseString, fileName, DefaultOpenMode<SparseString<TFibreValueString, TSpec> >::VALUE);
}

template <typename TFibreValueString, typename TSpec>
inline bool save(
    SparseString<TFibreValueString, TSpec> const & sparseString,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".val");
    if (!save(getFibre(sparseString, FibreValueString()), toCString(name), openMode))
    {
        return false;
    }
    name = fileName;    append(name, ".ind");   save(getFibre(sparseString, FibreIndicatorString()), toCString(name), openMode);
    return true;
}
// TODO(singer): setValue function

}
#endif // INDEX_FM_SPARSE_STRING_H_
