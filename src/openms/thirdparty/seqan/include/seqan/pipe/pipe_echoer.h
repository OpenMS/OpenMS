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

#ifndef SEQAN_HEADER_PIPE_ECHOER_H
#define SEQAN_HEADER_PIPE_ECHOER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    //////////////////////////////////////////////////////////////////////////////
	// some metaprogramming to unrool fixed-size loops
    struct EchoerFillWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.tmp.i2[I-1] = *(arg.in); ++(arg.in);
        }
    };
    
    struct EchoerClearWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
			arg.i2[I] = typename Value< typename Value<Arg, 2>::Type >::Type ();
        }
    };
    
    struct EchoerShiftWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.i2[I] = arg.i2[I-1];
        }
    };


    template < unsigned echoRepeats, bool omitFirst >
    struct Echoer;

    template < typename TInput, unsigned echoRepeats, bool omitFirst >
    struct Value< Pipe< TInput, Echoer< echoRepeats, omitFirst > > > {
        typedef Tuple<typename Value<TInput>::Type, echoRepeats>	EchoType;
        typedef Pair<typename Size<TInput>::Type, EchoType>			Type;
    };


/**
.Spec.Echoer:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs tuples of the $echoRepeats$ last elements of the input stream.
..signature:Pipe<TInput, Echoer<echoRepeats, omitFirst> >
..param.TInput:The type of the pipeline module this module reads from.
..param.echoRepeats:The tuple length.
...remarks:The tuples contain elements $in[i]in[i-1]...in[i-(echoRepeats-1)]$.
..param.omitFirst:Omit half filled tuples.
..param.omitFirst:If $true$, the output stream is $echoRepeats-1$ elements shorter than the input stream.
..param.omitFirst:If $false$, the lengths are identical and the tuple is filled with blanks (default constructed elements) for undefined entries.
..remarks:The output type is a @Class.Tuple@ of input elements and length $echoRepeats$ (i.e. $Tuple<Value<TInput>::Type, echoRepeats>$).
..remarks:The tuples are sequences of the form $in[i]in[i-1]in[i-2]..in[i-echoRepeats+1]$. For $omitFirst=false$ $i$ begins with 0 and for $omitFirst=true$ $i$ begins with $echoRepeats-1$.
..include:seqan/pipe.h
*/

    //////////////////////////////////////////////////////////////////////////////
    // echoer class
    template < typename TInput, unsigned echoRepeats, bool omitFirst >
    struct Pipe< TInput, Echoer<echoRepeats, omitFirst> >
    {
        typedef typename Value<Pipe>::Type TValue;

        TInput	&in;
        TValue	tmp;

        Pipe(TInput& _in):
            in(_in),
            tmp(0, typename Value<TValue, 2>::Type()) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			++in;
            if (eof(in)) return *this;
            LoopReverse<EchoerShiftWorker_, echoRepeats - 1>::run(this->tmp);
			++tmp.i1;
            tmp.i2[0] = *in;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned echoRepeats, bool omitFirst >
	inline bool control(Pipe< TInput, Echoer< echoRepeats, omitFirst > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.tmp.i1 = 0;
        Loop<EchoerClearWorker_, echoRepeats - 1>::run(me.tmp);
        if (!eof(me.in)) me.tmp.i2[0] = *me.in;
		return true;
	}
    
    template < typename TInput, unsigned echoRepeats >
    inline bool control(Pipe< TInput, Echoer< echoRepeats, true > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command) || size(me.in) < echoRepeats - 1) return false;
        me.tmp.i1 = 0;
        LoopReverse<EchoerFillWorker_, echoRepeats - 1>::run(me);
        if (!eof(me.in)) me.tmp.i2[0] = *me.in;
		return true;
    }

    template < typename TInput, unsigned echoRepeats >
    inline Size< Pipe< TInput, Echoer< echoRepeats, true > > >
    length(Pipe< TInput, Echoer< echoRepeats, true > > const &me) {
        return length(me.in) - (echoRepeats - 1);
    }

//}

}

#endif
