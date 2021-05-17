// ==========================================================================
//                            misc_accumulators.h
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
// Implementation of the SeqAn accumulator code.
//
// While Boost also has a versatile accumulator library, we create our own
// since it is much simpler and focused on the requirements in sequence
// analysis algorithms than then Boost one.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_MISC_MISC_ACCUMULATORS_H_
#define CORE_INCLUDE_SEQAN_MISC_MISC_ACCUMULATORS_H_

#include <seqan/sequence.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T, typename TTag> struct Result;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Average_;
typedef Tag<Average_> Average;

struct Sum_;
typedef Tag<Sum_> Sum;

struct Count_;
typedef Tag<Count_> Count;

/**
.Class.Accumulator
..cat:Miscellaneous
..summary:Accumulator base class.
..signature:Accumulator<TValue, TSpec>
..param.TValue:The type of the values to accumulate.
..param.TSpec:The specialization tag.
..remarks:Accumulators are for computing statistics on streams of values.
..remarks:Currently, this is only meant for accumulating integers.
..include:seqan/misc/misc_accumulators.h
*/

template <typename TValue, typename TSpec>
struct Accumulator;

/**
.Spec.Average Accumulator
..general:Class.Accumulator
..cat:Miscellaneous
..summary:Accumulator for computing averages.
..signature:Accumulator<TValue, TSpec>
..param.TValue:The type of the values to compute the average of.
..param.TSpec:The specialization tag.
..remarks:The average of an empty sequence is defined to be 0.
..example.text:This program shows how to use the Average Accumulator.
..example.code:
Accumulator<int, Average> acc;
push(acc, 1);
push(acc, 2);
push(acc, 3);
std::cout << "average: " << average(acc) << "\n"
          << "sum:     " << sum(acc) << "\n"
          << "count:   " << count(acc) << "\n";
..example.text:The output is then:
..example.code:
average: 2
sum:     6
count:   3
..include:seqan/misc/misc_accumulators.h
*/

template <typename TValue>
struct Accumulator<TValue, Average>
{
    typedef typename Result<Accumulator<TValue, Average>, Sum>::Type TSum_;

    TSum_ sum_;
    unsigned count_;

    Accumulator() : sum_(0), count_(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value                                       Average Accumulator
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Accumulator
///.Metafunction.Value.class:Class.Accumulator

template <typename TValue>
struct Value<Accumulator<TValue, Average> >
{
    typedef double Type;
};

template <typename TValue>
struct Value<Accumulator<TValue, Average> const > : Value<Accumulator<TValue, Average> >
{};

// ----------------------------------------------------------------------------
// Metafunction Result                                      Average Accumulator
// ----------------------------------------------------------------------------

// TODO(holtgrew): This could probably go to basic, as a part of an "AlgorithmState" concept?

/**
.Metafunction.Result
..cat:Miscellaneous
..summary:Return the result of a computation.
..signature:Result<T, TTag>::Type
..param.T:The type to query
..param.TTag:The type of the result to query.
...remarks:E.g. $Average$, $Sum$, $Count$.
..include:seqan/misc/misc_accumulators.h
 */

template <typename T, typename TTag = void>
struct Result;

// ----------------------------------------------------------------------------
// Metafunction Result                                      Average Accumulator
// ----------------------------------------------------------------------------

///.Metafunction.Result.param.T.type:Class.Accumulator

template <typename TValue>
struct Result<Accumulator<TValue, Average>, Average>
{
    typedef double Type;
};

template <typename TValue>
struct Result<Accumulator<TValue, Average> const, Average> : Result<Accumulator<TValue, Average>, Average>
{};

template <typename TValue>
struct Result<Accumulator<TValue, Average>, Count>
{
    typedef unsigned Type;
};

template <typename TValue>
struct Result<Accumulator<TValue, Average> const, Count> : Result<Accumulator<TValue, Average>, Count>
{};

template <typename TValue>
struct Result<Accumulator<TValue, Average>, Sum>
{
    typedef typename IfC<Is<IntegerConcept<TValue> >::VALUE,
                         __int64,
                         double>::Type Type;
};

template <typename TValue>
struct Result<Accumulator<TValue, Average> const, Sum> : Result<Accumulator<TValue, Average>, Sum>
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()                                                 Accumulator
// ----------------------------------------------------------------------------

///.Function.clear.param.object.type:Class.Accumulator

// ----------------------------------------------------------------------------
// Function push()                                                  Accumulator
// ----------------------------------------------------------------------------

/**
.Function.Accumulator#push
..cat:Miscellaneous
..summary:Adds a value to an accumulator.
..signature:push(acc, x)
..class:Class.Accumulator
..param.acc:The accumulator to add the value to.
...type:Class.Accumulator
..param.x:The value to push.
..returns:$void$
*/

// ----------------------------------------------------------------------------
// Function clear()                                         Average Accumulator
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
clear(Accumulator<TValue, Average> & accumulator)
{
    accumulator.sum_ = 0;
    accumulator.count_ = 0;
}

// ----------------------------------------------------------------------------
// Function push()                                          Average Accumulator
// ----------------------------------------------------------------------------

template <typename TValue, typename TValue2>
inline void
push(Accumulator<TValue, Average> & acc, TValue2 value)
{
    typedef typename Result<Accumulator<TValue, Average>, Sum>::Type TSum;
    acc.sum_ += static_cast<TSum>(value);
    acc.count_ += 1;
}

// ----------------------------------------------------------------------------
// Function average()                                       Average Accumulator
// ----------------------------------------------------------------------------

/**
.Function.Accumulator#average
..cat:Miscellaneous
..summary:Return average from an accumulator.
..signature:average(acc)
..class:Spec.Average Accumulator
..param.acc:The accumulator to return the average from.
...type:Spec.Average Accumulator
..returns:nolink:$double$
..include:seqan/misc/misc_accumulators.h
*/

template <typename TValue>
inline typename Result<Accumulator<TValue, Average>, Average>::Type
average(Accumulator<TValue, Average> const & acc)
{
    typedef typename Result<Accumulator<TValue, Average>, Average>::Type TResult;
    if (acc.count_ == 0u)
        return 0;
    return static_cast<TResult>(acc.sum_ / static_cast<TResult>(acc.count_));
}

// ----------------------------------------------------------------------------
// Function sum()                                           Average Accumulator
// ----------------------------------------------------------------------------

/**
.Function.Accumulator#sum
..cat:Miscellaneous
..summary:Return sum from an accumulator.
..signature:average(acc)
..class:Spec.Average Accumulator
..param.acc:The accumulator to return the sum from.
...type:Spec.Average Accumulator
..returns:@Metafunction.Value@$<typeof(acc)>::Type$
..include:seqan/misc/misc_accumulators.h
*/

template <typename TValue>
inline typename Result<Accumulator<TValue, Average>, Sum>::Type
sum(Accumulator<TValue, Average> const & acc)
{
    return acc.sum_;
}

// ----------------------------------------------------------------------------
// Function count()                                         Average Accumulator
// ----------------------------------------------------------------------------

/**
.Function.Accumulator#count
..cat:Miscellaneous
..summary:Return sum from an accumulator.
..signature:average(acc)
..class:Spec.Average Accumulator
..param.acc:The accumulator to return the sum from.
...type:Spec.Average Accumulator
..returns:@Metafunction.Value@$<typeof(acc)>::Type$
..include:seqan/misc/misc_accumulators.h
*/

template <typename TValue>
inline typename Result<Accumulator<TValue, Average>, Count>::Type
count(Accumulator<TValue, Average> const & acc)
{
    return acc.count_;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_MISC_MISC_ACCUMULATORS_H_
