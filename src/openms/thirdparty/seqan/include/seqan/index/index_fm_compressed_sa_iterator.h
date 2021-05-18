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

#ifndef INDEX_FM_COMPRESSED_SA_ITERATOR_H_
#define INDEX_FM_COMPRESSED_SA_ITERATOR_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

struct FibreSparseString_;
typedef Tag<FibreSparseString_> const FibreSparseString;

template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA;

// ==========================================================================
// Metafunctions
// ==========================================================================

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>
{
    typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec> const, PositionIterator> Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>
{
    typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec>, PositionIterator> Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>:
    Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>{};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>:
    Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>{};

// ==========================================================================
// Functions
// ==========================================================================
///.Function.begin.param.object.type:Class.CompressedSA
template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>::Type
begin(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, Standard const & /* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>::Type(compressedSA, 0);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>::Type
begin(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, Standard const & /* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>::Type(compressedSA, 0);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>::Type
begin(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, Rooted const & /* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>::Type(compressedSA, 0);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>::Type
begin(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, Rooted const & /* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>::Type(compressedSA, 0);
}

// ==========================================================================
///.Function.end.param.object.type:Class.CompressedSA
template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>::Type
end(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, Rooted const & /* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>::Type(compressedSA, length(compressedSA));
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>::Type
end(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, Rooted/* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>::Type(compressedSA, length(compressedSA));
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>::Type
end(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, Standard const & /* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>::Type(compressedSA, length(compressedSA));
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>::Type
end(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, Standard/* dummy */)
{
    return typename Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>::Type(compressedSA, length(compressedSA));
}

}
#endif // INDEX_FM_COMPRESSED_SA_ITERATOR_H_
