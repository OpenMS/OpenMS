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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_PIPE_FILTER_H
#define SEQAN_HEADER_PIPE_FILTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{
    
    template <typename TValue, typename TResult = typename Value<TValue, 1>::Type>
    struct filterI1 : public ::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue & x) const
        {
            return x.i1;
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 2>::Type>
    struct filterI2 : public ::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue & x) const
        {
            return x.i2;
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 3>::Type>
    struct filterI3 : public ::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue & x) const
        {
            return x.i3;
        }
    };


    template < typename TFunctor >
    struct Filter;

	template < typename TInput, typename TFunctor >
    struct Value< Pipe< TInput, Filter<TFunctor> > >
    {
		typedef typename TFunctor::result_type Type;
	};


/**
.Spec.Filter:
..cat:Pipelining
..general:Class.Pipe
..summary:Applies a specific function to the input stream.
..signature:Pipe<TInput, Filter<TFunctor> >
..param.TInput:The type of the pipeline module this module reads from.
..param.TFunctor:A unary function (see STL's $unary_function$).
...remarks:The argument type of $TFunctor$ must be $VALUE<TInput>::Type$.
..remarks: The output type of this pipe is the result type of $TFunctor$.
..include:seqan/pipe.h
*/

	//////////////////////////////////////////////////////////////////////////////
    // filter class
    template <typename TInput, typename TFunctor >
    struct Pipe< TInput, Filter<TFunctor> >
    {
		TInput      &in;
        TFunctor    F;
        
/**
.Memfunc.Filter#Pipe:
..class:Spec.Filter
..summary:Constructor
..signature:Pipe<TInput, Filter<TFunctor> > (in)
..signature:Pipe<TInput, Filter<TFunctor> > (in, func)
..param.in:Reference to an input pipe.
..param.func:A TFunctor object (copy constructor).
*/
        Pipe(TInput& _in):
            in(_in) {}
        
        Pipe(TInput& _in, const TFunctor& F_) :
            in(_in),
            F(F_) {}
        
        inline typename Value<Pipe>::Type const operator*() const
        {
            return F(*in);
        }

        Pipe & operator++()
        {
            ++in;
            return *this;
        }
                
    };
    
//}

}

#endif
