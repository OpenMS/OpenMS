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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Various useful bit-twiddling routines, mostly taken from the website
// http://www-graphics.stanford.edu/~seander/bithacks.html
// ==========================================================================

#ifndef SEQAN_MISC_MISC_BIT_TWIDDLING_H_
#define SEQAN_MISC_MISC_BIT_TWIDDLING_H_

#ifdef PLATFORM_WINDOWS_VS
// Make intrinsics visible.  It appears that this is not necessary with VS 10
// any more, for VS 9, it must be included.
#include <intrin.h>
#endif  // #ifdef PLATFORM_WINDOWS_VS

// TODO(holtgrew): Test this!

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setBitTo()
// ----------------------------------------------------------------------------

/**
.Function.setBitTo
..cat:Bit Twiddling
..summary:Set the bit with the given index to the given value.
..signature:setBitTo(word, index, value)
..param.word:The number.
..param.index:The index of the bit in the word.
...type:nolink:$unsigned$
..param.value:The value to set the bit to.
...type:nolink:$bool
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBit
..see:Function.clearBit
..see:Function.clearAllBits
..see:Function.isBitSet
 */

template <typename TWord, typename TPos>
inline void
setBitTo(TWord & word, TPos index, bool value)
{
    // See http://www-graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
    word = (word & ~(1u << index)) | (-value & (1u << index));
}

// ----------------------------------------------------------------------------
// Function setBit()
// ----------------------------------------------------------------------------

/**
.Function.setBit
..cat:Bit Twiddling
..summary:Set the bit with the given index to 1.
..signature:setBit(word, index)
..param.word:The number.
..param.index:The index of the bit in the word.
...type:nolink:$unsigned$
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.clearBit
..see:Function.clearAllBits
..see:Function.isBitSet
 */

template <typename TWord, typename TPos>
inline void
setBit(TWord & word, TPos index)
{
    word |= (1u << index);
}

// ----------------------------------------------------------------------------
// Function clearBit()
// ----------------------------------------------------------------------------

/**
.Function.clearBit
..cat:Bit Twiddling
..summary:Set the bit with the given index to 0.
..signature:clearBit(word, index)
..param.word:The number.
..param.index:The index of the bit in the word.
...type:nolink:$unsigned$
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.setBit
..see:Function.clearAllBits
..see:Function.isBitSet
 */

template <typename TWord, typename TPos>
inline void
clearBit(TWord & word, TPos index)
{
    word &= ~(1u << index);
}

// ----------------------------------------------------------------------------
// Function clearAllBits()
// ----------------------------------------------------------------------------

/**
.Function.clearAllBits
..cat:Bit Twiddling
..summary:Set all bits to 0.
..signature:clearAllBits(word)
..param.word:The number.
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.setBit
..see:Function.clearBit
..see:Function.isBitSet
 */

template <typename TWord>
inline
void
clearBits(TWord & word)
{
    word = 0;
}

// ----------------------------------------------------------------------------
// Function isBitSet()
// ----------------------------------------------------------------------------

/**
.Function.isBitSet
..cat:Bit Twiddling
..summary:Returns whether the bit with the given index is set to 1.
..signature:isBitSet(word, index)
..param.word:The number.
..param.index:The index.
...type:nolink:$unsigned$
..returns:$true$ if the bit was set and $false$ otherwise.
...type:nolink:$bool$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.setBit
..see:Function.clearBit
..see:Function.clearAllBits
 */

template <typename TWord, typename TIndex>
inline
bool
isBitSet(TWord const & word, TIndex index)
{
    typedef typename MakeUnsigned<TWord>::Type TUnsignedWord;
    return (word & (TUnsignedWord(1) << index)) != static_cast<TWord>(0);
}

// ----------------------------------------------------------------------------
// Function popCount()
// ----------------------------------------------------------------------------

/**
.Function.popCount
..cat:Bit Twiddling
..summary:Returns number of set bits in an integer.
..signature:popCount(word)
..param.word:The number.
..returns:The number of set bits (1s) in an integer.
...type:nolink:$unsigned$
..include:seqan/misc/misc_bit_twiddling.h
 */

// Implementing this platform-independent is tricky.  There are two points to platform independentness. First, the
// choice of compiler and second the used CPU.  Currently, we do not perform any checks for the CPU and assume that
// the Intel intrinsic POPCNT is available.  The function is implemented to work on the supported compilers GCC/MINGW,
// CLANG (which has the same interface as GCC here) and Visual C++.
//
// GCC, MINGW and CLANG provide the intrinsics __builtin_popcount, __builtin_popcountl, and __builtin_popcountll for
// the types unsigned, unsigned long, and unsigned long long.  Starting with version 2008, Visual C++ provides the
// intrinsics __popcnt16, __popcnt, and __popcnt64 for 16, 32, and 64 bit words.
//
// The functions below are implemented as follows.  _popCountImplGeneric() is used if there are no intrinsics provided
// by the compiler (the case for Visual C++ 2008).  Otherwise, we define different overloads of the function
// _popCountImpl() that are given the length of the word as a template argument.  If necessary, we copy the word in a
// variable of next largest size and call the best suited builtin on this copy.

// Generic implementation of counting bits.  Taken from http://graphics.stanford.edu/~seander/bithacks.html
//
// Brian Kernighan's method goes through as many iterations as there are set bits. So if we have a 32-bit word with
// only the high bit set, then it will only go once through the loop.
//
// Published in 1988, the C Programming Language 2nd Ed. (by Brian W. Kernighan and Dennis M. Ritchie) mentions this
// in exercise 2-9. On April 19, 2006 Don Knuth pointed out to me that this method "was first published by Peter
// Wegner in CACM 3 (1960), 322. (Also discovered independently by Derrick Lehmer and published in 1964 in a book
// edited by Beckenbach.)"

template <typename TWord>
inline unsigned
_popCountImplGeneric(TWord const & word)  // Note that word is copied!
{
    typename MakeUnsigned<TWord>::Type x = word;
	unsigned int c = 0;  // c accumulates the total bits set in v
	for (c = 0; x; c++)
		x &= x - 1;  // clear the least significant bit set
	return c;
}

// This parametrized tag is used for selecting a _popCountImpl() implementation.

template <unsigned int NUM_BITS>
struct WordSize_ {};

// The compiler-dependent implementations of _popCountImpl() follow.

#if defined(_MSC_VER) && (_MSC_VER <= 1400)  // MSVC <= 2005, no intrinsic.

template <typename TWord, unsigned NUM_BITS>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<NUM_BITS> const & /*tag*/)
{
    return _popCountImplGeneric(word);
}

#endif  // #if defined(_MSC_VER) && (_MSC_VER <= 1400)  // MSVC <= 2005, no intrinsic.

#if defined(_MSC_VER) && (_MSC_VER > 1400)  // MSVC >= 2008, has intrinsic

#if defined(_WIN64)

// 64-bit Windows, 64 bit intrinsic available

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<64> const & /*tag*/)
{
    return __popcnt64(static_cast<__uint64>(word));
}

#else  // #if defined(_WIN64)

// 32-bit Windows, 64 bit intrinsic not available

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<64> const & /*tag*/)
{
	return __popcnt(static_cast<__uint32>(word)) + __popcnt(static_cast<__uint32>(word >> 32));
}

#endif  // #if defined(_WIN64)

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<32> const & /*tag*/)
{
    return __popcnt(static_cast<__uint32>(word));
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<16> const & /*tag*/)
{
    return __popcnt16(static_cast<__uint16>(word));
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<8> const & /*tag*/)
{
    return _popCountImpl(static_cast<const __uint16>(word), WordSize_<16>());
}

#endif  // #if defined(_MSC_VER) && (_MSC_VER <= 1400)

#if !defined(_MSC_VER)  // GCC or CLANG

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<64> const & /*tag*/)
{
    return __builtin_popcountll(static_cast<unsigned long long>(word));
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<32> const & /*tag*/)
{
    return __builtin_popcount(static_cast<unsigned int>(word));
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<16> const & /*tag*/)
{
    return _popCountImpl(static_cast<__uint32>(word), WordSize_<32>());
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord const & word, WordSize_<8> const & /*tag*/)
{
    return _popCountImpl(static_cast<__uint32>(word), WordSize_<32>());
}

#endif    // GCC or CLANG

template <typename TWord>
inline unsigned
popCount(TWord word)
{
    return _popCountImpl(word, WordSize_<BitsPerValue<TWord>::VALUE>());
}

}  // namespace seqan

#endif // #ifndef SEQAN_MISC_MISC_BIT_TWIDDLING_H_
