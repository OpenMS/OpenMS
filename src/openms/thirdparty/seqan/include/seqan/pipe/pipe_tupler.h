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

#ifndef SEQAN_HEADER_PIPE_TUPLER_H
#define SEQAN_HEADER_PIPE_TUPLER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

//////////////////////////////////////////////////////////////////////////////

    template < unsigned tupleLen, bool omitLast = false, typename TPack = void >
    struct Tupler;

    template < typename TInput, unsigned tupleLen, bool omitLast, typename TPack >
    struct Value< Pipe< TInput, Tupler< tupleLen, omitLast, TPack > > >
    {
        typedef Tuple<typename Value<TInput>::Type, tupleLen, TPack>	TTuple;
        typedef Pair<typename Size<TInput>::Type, TTuple, Pack>         Type;
    };

//////////////////////////////////////////////////////////////////////////////

    template < 
		typename TInput, 
		unsigned tupleLen, 
		bool omitLast, 
		typename TPack,
		typename TPair, 
		typename TLimitsString >
    struct Value< Pipe< TInput, Multi< Tupler< tupleLen, omitLast, TPack >, TPair, TLimitsString > > >
    {
        typedef Tuple<typename Value<TInput>::Type, tupleLen, TPack>	TTuple;
        typedef Pair<TPair, TTuple, Pack>								Type;
    };

//////////////////////////////////////////////////////////////////////////////


	// output only fully filled tuples
	template < typename TTupler >
	struct TuplerLastTuples_ {
		enum { VALUE = 1 };
	};

	// output tupleLen-1 half filled tuples at the end
    template < typename TInput, unsigned tupleLen, typename TPack >
	struct TuplerLastTuples_< Pipe< TInput, Tupler<tupleLen, false, TPack> > > {
		enum { VALUE = tupleLen };
	};

    struct ShiftLeftWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.i2[I-1] = arg.i2[I];
        }
    };

/**
.Spec.Tupler:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs tuples of the $tupleLen$ consecutive elements of the input stream.
..signature:Pipe<TInput, Tupler<tupleLen, omitLast> >
..param.TInput:The type of the pipeline module this module reads from.
..param.tupleLen:The tuple length.
...remarks:The tuples contain elements $in[i]in[i+1]...in[i+(tupleLen-1)]$.
..param.omitLast:Omit half filled tuples.
..param.omitLast:If $true$, the output stream is $tupleLen-1$ elements shorter than the input stream.
..param.omitLast:If $false$, the lengths are identical and the last tuples are filled with blanks (default constructed elements) for undefined entries.
..remarks:The output type is a @Class.Tuple@ of input elements and length $tupleLen$ (i.e. $Tuple<Value<TInput>::Type, tupleLen>$).
..remarks:The tuples are sequences of the form $in[i]in[i-1]in[i-2]..in[i-tupleLen+1]$. For $omitLast=false$ $i$ begins with 0 and for $omitLast=true$ $i$ begins with $tupleLen-1$.
..include:seqan/pipe.h
*/

    //////////////////////////////////////////////////////////////////////////////
    // tupler class
    template < typename TInput, unsigned tupleLen, bool omitLast, typename TPack >
    struct Pipe< TInput, Tupler<tupleLen, omitLast, TPack> >
    {
		typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		TInput                      &in;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	lastTuples;
        
        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
            if (eof(in)) --lastTuples;
            Loop<ShiftLeftWorker_, tupleLen - 1>::run(this->tmp);
			++tmp.i1;
			if (lastTuples < TuplerLastTuples_<Pipe>::VALUE)
	            tmp.i2[tupleLen - 1] = TValue();
			else {
				tmp.i2[tupleLen - 1] = *in;
				++in;
			}
            return *this;
        }

        inline void fill() {
            unsigned i;
            for(i = 0; i < tupleLen && !eof(in); ++i, ++in)
                tmp.i2.i[i] = *in;
			if (TuplerLastTuples_<Pipe>::VALUE > tupleLen - i)
				lastTuples = TuplerLastTuples_<Pipe>::VALUE - (tupleLen - i);
			else
				lastTuples = 0;
            for(; i < tupleLen; ++i)
                tmp.i2.i[i] = TValue();
            tmp.i1 = 0;
        }
	};

//____________________________________________________________________________


	template < typename TInput, unsigned tupleLen, bool omitLast >
    struct Pipe< TInput, Tupler<tupleLen, omitLast, BitPacked<> > >
    {
        TInput                      &in;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	lastTuples;
        
        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const
        {
            return tmp;
        }

        inline Pipe& operator++()
        {
            if (eof(in)) --lastTuples;
			tmp.i2 <<= 1;
			++tmp.i1;
			if (lastTuples == TuplerLastTuples_<Pipe>::VALUE)
            {
				tmp.i2 |= *in;
				++in;
			}
            return *this;
        }

        inline void fill()
        {
            unsigned i;
            clear(tmp.i2);
			for(i = 0; i < tupleLen && !eof(in); ++i, ++in)
            {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
			}
			if (TuplerLastTuples_<Pipe>::VALUE > tupleLen - i)
				lastTuples = TuplerLastTuples_<Pipe>::VALUE - (tupleLen - i);
			else
				lastTuples = 0;
            tmp.i2 <<= (tupleLen - i);
            tmp.i1 = 0;
        }
	};


    //////////////////////////////////////////////////////////////////////////////
    // tupler class for multiple sequences
    template < 
		typename TInput, 
		unsigned tupleLen, 
		bool omitLast, 
		typename TPack, 
		typename TPair, 
		typename TLimitsString >
    struct Pipe< TInput, Multi<Tupler<tupleLen, omitLast, TPack>, TPair, TLimitsString> >
    {
		typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		typedef PairIncrementer_<TPair, TLimitsString>	Incrementer;

		TInput                      &in;
        Incrementer					localPos;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	seqLength, lastTuples;

		TLimitsString const &limits;

        template <typename TLimitsString_>
        Pipe(TInput& _in, TLimitsString_ &_limits):  // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
            in(_in),
			limits(_limits) {}

        inline typename Value<Pipe>::Type const & operator*() const
        {
            return tmp;
        }

        inline Pipe& operator++()
        {
			// process next sequence
			if (eos())
				if (--lastTuples == 0)
                {
					assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);
					fill();
					return *this;
				}

			// shift left 1 character
            Loop<ShiftLeftWorker_, tupleLen - 1>::run(this->tmp);
			assignValueI2(tmp.i1, getValueI2(tmp.i1) + 1);

			if (lastTuples < TuplerLastTuples_<Pipe>::VALUE)
            {
	            tmp.i2[tupleLen - 1] = TValue();
			} else
            {
				tmp.i2[tupleLen - 1] = *in;
				++localPos;
				++in;
			}
            return *this;
        }

        inline void fill()
        {
			do {
				unsigned i = 0;
				if (!eof(in))
					do {
						tmp.i2.i[i] = *in;
						++in;
						++i;
						++localPos;
					} while ((i < tupleLen) && !eos());
				lastTuples = TuplerLastTuples_<Pipe>::VALUE;

				// fill up with null chars
				for(; i < tupleLen; ++i)
					tmp.i2.i[i] = TValue();
				
				// eventually, reduce the number of half-filled tuples
				if (lastTuples <= tupleLen - i)
					lastTuples = 0;
				else
					lastTuples -= tupleLen - i;

				if (lastTuples == 0)
					assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);

			} while ((lastTuples == 0) && !eof(in));

			assignValueI2(tmp.i1, 0);
		}

		inline bool eos()
        {
			return (getValueI1(localPos) > 0) && (getValueI2(localPos) == 0);
		}
	};

//____________________________________________________________________________


	template < 
		typename TInput, 
		unsigned tupleLen, 
		bool omitLast, 
		typename TPair, 
		typename TLimitsString >
    struct Pipe< TInput, Multi<Tupler<tupleLen, omitLast, BitPacked<> >, TPair, TLimitsString> >
    {
		typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		typedef PairIncrementer_<TPair, TLimitsString>	Incrementer;

		TInput                      &in;
        Incrementer					localPos;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	seqLength, lastTuples;

		TLimitsString const &limits;
        
        template <typename TLimitsString_>
        Pipe(TInput& _in, TLimitsString_ &_limits):  // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
            in(_in),
			limits(_limits) {}

        inline typename Value<Pipe>::Type const & operator*() const
        {
            return tmp;
        }

        inline Pipe& operator++() {
			// process next sequence
			if (eos())
				if (--lastTuples == 0) {
					assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);
					fill();
					return *this;
				}

			// shift left 1 character
			tmp.i2 <<= 1;
			assignValueI2(tmp.i1, getValueI2(tmp.i1) + 1);
			if (lastTuples == TuplerLastTuples_<Pipe>::VALUE) {
				tmp.i2 |= *in;
				++localPos;
				++in;
			}
            return *this;
        }

        inline void fill()
        {
			do
            {
				unsigned i = 0;
				if (!eof(in))
					do
                    {
						tmp.i2 <<= 1;
						tmp.i2 |= *in;
						++in;
						++i;
						++localPos;
					} while ((i < tupleLen) && !eos());
				lastTuples = TuplerLastTuples_<Pipe>::VALUE;

				// fill up with null chars
	            tmp.i2 <<= (tupleLen - i);
				
				// eventually, reduce the number of half-filled tuples
				if (lastTuples <= tupleLen - i)
					lastTuples = 0;
				else
					lastTuples -= tupleLen - i;

				if (lastTuples == 0)
					assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);

			} while ((lastTuples == 0) && !eof(in));

			assignValueI2(tmp.i1, 0);
        }

		inline bool eos()
        {
			return (getValueI1(value(localPos)) > 0) && (getValueI2(value(localPos)) == 0);
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned tupleLen, bool omitLast, typename TPack >
	inline bool 
	control(
		Pipe< TInput, Tupler< tupleLen, omitLast, TPack > > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;
		me.fill();
		return true;
	}
    
    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TPack,
		typename TPair, 
		typename TLimitsString >
	inline bool 
	control(
		Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TPack >, TPair, TLimitsString> > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;
		setHost(me.localPos, me.limits);
		assignValueI1(me.tmp.i1, 0);
		me.fill();
		return true;
	}
    
    template < typename TInput, unsigned tupleLen, bool omitLast, typename TPack >
	inline bool 
	control(
		Pipe< TInput, Tupler< tupleLen, omitLast, TPack > > &me, 
		ControlEof const &)
	{
		return me.lastTuples == 0;
    }

    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TPack,
		typename TPair, 
		typename TLimitsString >
	inline bool 
	control(
		Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TPack >, TPair, TLimitsString> > &me, 
		ControlEof const &) 
	{
		return me.lastTuples == 0;
	}

    template < typename TInput, unsigned tupleLen, bool omitLast, typename TPack >
	inline bool 
	control(
		Pipe< TInput, Tupler< tupleLen, omitLast, TPack > > &me, 
		ControlEos const &) 
	{
		return control(me, ControlEof());
	}

    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TPack,
		typename TPair, 
		typename TLimitsString >
	inline bool 
	control(
		Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TPack >, TPair, TLimitsString> > &me, 
		ControlEos const &) 
	{
		return (getValueI1(me.tmp.i1) > 0) && (getValueI2(me.tmp.i1) == 0);
	}

    template < typename TInput, unsigned tupleLen, bool omitLast, typename TPack >
    inline typename Size< Pipe< TInput, Tupler< tupleLen, omitLast, TPack > > >::Type
    length(Pipe< TInput, Tupler< tupleLen, omitLast, TPack > > const &me) 
	{
		typedef Pipe< TInput, Tupler< tupleLen, omitLast, TPack > >	TPipe;
		if (length(me.in) >= (tupleLen - TuplerLastTuples_<TPipe>::VALUE))
			return length(me.in) - (tupleLen - TuplerLastTuples_<TPipe>::VALUE);
		else
			return 0;
    }

    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TPack,
		typename TPair, 
		typename TLimitsString >
    inline typename Size< Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TPack >, TPair, TLimitsString> > >::Type
    length(Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TPack >, TPair, TLimitsString> > const &me)
	{
		typedef Pipe< TInput, Tupler< tupleLen, omitLast, TPack > >	TPipe;
		unsigned seqs = countSequences(me);
		
		if (length(me.in) >= seqs * (tupleLen - TuplerLastTuples_<TPipe>::VALUE))
			return length(me.in) - seqs * (tupleLen - TuplerLastTuples_<TPipe>::VALUE);
		else
			return 0;
    }

    template < typename TInput, unsigned tupleLen, bool omitLast, typename TPack >
    inline unsigned
    countSequences(Pipe< TInput, Tupler< tupleLen, omitLast, TPack > > const &)
    {
		return 1;
	}

    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TPack,
		typename TPair, 
		typename TLimitsString >
    inline unsigned
	countSequences(Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TPack >, TPair, TLimitsString> > const &me)
    {
		return length(me.limits) - 1;
	}

//}

}

#endif
