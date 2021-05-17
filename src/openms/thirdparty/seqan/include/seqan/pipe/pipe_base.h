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

#ifndef SEQAN_HEADER_PIPE_BASE_H
#define SEQAN_HEADER_PIPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

	// shortcuts to ease pipeline construction
    #define TypeOf_(TObject)  typename Value<TObject>::Type
    #define TSizeOf_(TObject) typename Size<TObject>::Type

/**
.Class.Pipe:
..cat:Pipelining
..summary:Pipes are pop-passive pipeline modules.
..signature:Pipe<TInput, TSpec>
..param.TInput:The type of the pipeline module this module reads from.
...remarks:Use @Class.Bundle2@, @Class.Bundle3@, etc. to read from more than one module.
..param.TSpec:The specializing type.
..remarks:Use @Metafunction.Value@ to get the output type of a given Pipe (returns $Value<TInput>::Type$ by default).
..remarks:Use @Metafunction.Size@ to get the size type of a given Pipe (returns $Size<TInput>::Type$ by default).
.Memfunc.Pipe#Pipe:
..class:Class.Pipe
..summary:Constructor
..signature:Pipe<TInput, TSpec> (in)
..param.in:Reference to an input pipe.
..include:seqan/pipe.h
*/

    template < typename TInput, typename TSpec >
    struct Pipe {
        TInput &in;
        Pipe(TInput &_in): in(_in) {}
    };

	// base class for multiple sequence algorithms
	// to hold extra information about desired position type (TPair)
	// and the type storing absolute sequence offsets (TLimitsString)
    template <typename TSpec, typename TPair, typename TLimitsString>
	struct Multi;

/**
.Class.Bundle2:
..cat:Aggregates
..summary:Stores references to two arbitrary objects.
..signature:Bundle2<TInput1, TInput2>
..param.TInput1:The type of the first object.
..param.TInput2:The type of the second object.
..remarks:Primarily used as an adaptor for pipes with two sources.
.Memvar.Bundle2#in1:
..class:Class.Bundle2
..summary:TInput1 reference
.Memvar.Bundle2#in2:
..class:Class.Bundle2
..summary:TInput2 reference
..include:seqan/pipe.h
*/
    // pipe input adapter 2->1 pipe
    template < typename TInput1, typename TInput2 >
    struct Bundle2 {
        typedef TInput1 Input1;
        typedef TInput2 Input2;
        TInput1 &in1;
        TInput2 &in2;
        Bundle2(TInput1 &_in1, TInput2 &_in2): in1(_in1),in2(_in2) {}
    };

/**
.Function.bundle2
..class:Class.Bundle2
..cat:Pipelining
..summary:Returns a bundle of two objects.
..signature:bundle2(in1, in2)
..param.in1:First object.
..param.in2:Second object.
..returns:A @Class.Bundle2@ with references to $in1$ and $in2$.
..see:Class.Bundle2
..include:seqan/pipe.h
*/
	template < typename TInput1, typename TInput2 >
	inline Bundle2< TInput1, TInput2 >
	bundle2(TInput1 &_in1, TInput2 &_in2) {
		return Bundle2< TInput1, TInput2 >(_in1, _in2);
	}

/**
.Class.Bundle3:
..cat:Aggregates
..summary:Stores references to three arbitrary objects.
..signature:Bundle3<TInput1, TInput2, TInput3>
..param.TInput1:The type of the first object.
..param.TInput2:The type of the second object.
..param.TInput3:The type of the third object.
..remarks:Primarily used as an adaptor for pipes with three sources.
.Memvar.Bundle3#in1:
..class:Class.Bundle3
..summary:TInput1 reference
.Memvar.Bundle3#in2:
..class:Class.Bundle3
..summary:TInput2 reference
.Memvar.Bundle3#in3:
..class:Class.Bundle3
..summary:TInput3 reference
..include:seqan/pipe.h
*/
    // pipe input adapter 3->1 pipe
    template < typename TInput1, typename TInput2, typename TInput3 >
    struct Bundle3 {
        typedef TInput1 Input1;
        typedef TInput2 Input2;
        typedef TInput3 Input3;
        TInput1 &in1;
        TInput2 &in2;
        TInput3 &in3;
        Bundle3(TInput1 &_in1, TInput2 &_in2, TInput3 &_in3): in1(_in1),in2(_in2),in3(_in3) {}
    };

/**
.Function.bundle3
..class:Class.Bundle3
..cat:Pipelining
..summary:Returns a bundle of three objects.
..signature:bundle3(in1, in2, in3)
..param.in1:First object.
..param.in2:Second object.
..param.in3:Third object.
..returns:A @Class.Bundle3@ with references to $in1$, $in2$, and $in3$.
..see:Class.Bundle3
..include:seqan/pipe.h
*/
	template < typename TInput1, typename TInput2, typename TInput3 >
	inline Bundle3< TInput1, TInput2, TInput3 >
	bundle3(TInput1 &_in1, TInput2 &_in2, TInput3 &_in3) {
		return Bundle3< TInput1, TInput2, TInput3 >(_in1, _in2, _in3);
	}

/**
.Class.Bundle5:
..cat:Aggregates
..summary:Stores references to five arbitrary objects.
..signature:Bundle5<TInput1, TInput2, TInput3, TInput4, TInput5>
..param.TInput1:The type of the first object.
..param.TInput2:The type of the second object.
..param.TInput3:The type of the third object.
..param.TInput4:The type of the fourth object.
..param.TInput5:The type of the fifth object.
..remarks:Primarily used as an adaptor for pipes with five sources.
.Memvar.Bundle5#in1:
..class:Class.Bundle5
..summary:TInput1 reference
.Memvar.Bundle5#in2:
..class:Class.Bundle5
..summary:TInput2 reference
.Memvar.Bundle5#in3:
..class:Class.Bundle5
..summary:TInput3 reference
.Memvar.Bundle5#in4:
..class:Class.Bundle5
..summary:TInput4 reference
.Memvar.Bundle5#in5:
..class:Class.Bundle5
..summary:TInput5 reference
..include:seqan/pipe.h
*/
    // pipe input adapter 5->1 pipe
    template < typename TIn1, typename TIn2, typename TIn3, typename TIn4, typename TIn5 >
    struct Bundle5 {
        TIn1 &in1; TIn2 &in2;
        TIn3 &in3; TIn4 &in4;
        TIn5 &in5;
        Bundle5(TIn1& _in1, TIn2& _in2,
                TIn3& _in3, TIn4& _in4,
                TIn5& _in5):    in1(_in1),in2(_in2),
                                in3(_in3),in4(_in4),
                                in5(_in5) {}
    };

/**
.Function.bundle5
..class:Class.Bundle5
..cat:Pipelining
..summary:Returns a bundle of five objects.
..signature:bundle5(in1, in2, in3, in4, in5)
..param.in1:First object.
..param.in2:Second object.
..param.in3:Third object.
..param.in4:Fourth object.
..param.in5:Fifth object.
..returns:A @Class.Bundle5@ with references to $in1$, $in2$, $in3$, $in4$, and $in5$.
..see:Class.Bundle5
..include:seqan/pipe.h
*/
    template < typename TIn1, typename TIn2, typename TIn3, typename TIn4, typename TIn5 >
	inline Bundle5< TIn1, TIn2, TIn3, TIn4, TIn5 >
	bundle5(TIn1 &_in1, TIn2 &_in2, TIn3 &_in3, TIn4 &_in4, TIn5 &_in5) {
		return Bundle5< TIn1, TIn2, TIn3, TIn4, TIn5 >(_in1, _in2, _in3, _in4, _in5);
	}

    template < typename TValue, typename TSize >
    struct AbstractSource {};

    template < typename TValue, typename TSize >
    struct Value< Pipe<void, AbstractSource<TValue, TSize> > > {
        typedef TValue Type;
    };

    template < typename TValue, typename TSize >
    struct Size< Pipe<void, AbstractSource<TValue, TSize> > > {
        typedef TSize Type;
    };



    
    template < typename TInput, typename TSpec >
    struct Value< Pipe<TInput, TSpec> > {
        typedef typename Value<TInput>::Type Type;
    };

    template < typename TInput, typename TSpec >
    struct Size< Pipe<TInput, TSpec> > {
        typedef typename Size<TInput>::Type Type;
    };

    template < typename TInput1, typename TInput2 >
    struct Size< Bundle2< TInput1, TInput2 > > {
        typedef typename Size<TInput1>::Type Type;
    };

    template < typename TInput1, typename TInput2, typename TInput3 >
    struct Size< Bundle3< TInput1, TInput2, TInput3 > > {
        typedef typename Size<TInput1>::Type Type;
    };

    template < typename TIn1, typename TIn2, typename TIn3, typename TIn4, typename TIn5 >
    struct Size< Bundle5< TIn1, TIn2, TIn3, TIn4, TIn5 > > {
        typedef typename Size<TIn1>::Type Type;
    };

    template < typename TInput, typename TSpec >
    struct Position< Pipe<TInput, TSpec> > {
        typedef typename Size<Pipe<TInput, TSpec> >::Type Type;
    };

    template < typename TInput, typename TSpec >
    struct Difference< Pipe<TInput, TSpec> > {
		typedef typename MakeSigned_<typename Size<Pipe<TInput, TSpec> >::Type>::Type Type;
    };
/*
    template < typename TInput, typename TSpec >
	struct Iterator< Pipe<TInput, TSpec> >;

    template < typename TInput, typename TSpec >
	struct Iterator< Pipe<TInput, TSpec> const >:
		Iterator< Pipe<TInput, TSpec> > {};
*/

	template <typename T>
	struct Source;

	template <typename TInput, typename TSpec>
	struct Source<Pipe<TInput, TSpec> >
	{
		typedef TInput Type;
	};

	template < typename TInput, typename TSpec >
    inline TInput const &
    source(Pipe<TInput, TSpec> const &me) {
SEQAN_CHECKPOINT
        return me.in;
    }

	template < typename TInput, typename TSpec >
    inline TInput &
    source(Pipe<TInput, TSpec> &me) {
SEQAN_CHECKPOINT
        return me.in;
    }


///.Function.length.param.object.type:Class.Pipe
///.Function.length.class:Class.Pipe

	template < typename TInput, typename TSpec >
    inline typename Size< Pipe<TInput, TSpec> >::Type
    length(Pipe<TInput, TSpec> const &me) {
SEQAN_CHECKPOINT
        return length(me.in);
    }

    template < typename TInput1, typename TInput2 >
    inline typename Size< Bundle2<TInput1, TInput2> >::Type
    length(Bundle2<TInput1, TInput2> const &me) {
SEQAN_CHECKPOINT
        return length(me.in1);
    }

    template < typename TInput1, typename TInput2, typename TInput3 >
    inline typename Size< Bundle3<TInput1, TInput2, TInput3> >::Type
    length(Bundle3<TInput1, TInput2, TInput3> const &me) {
SEQAN_CHECKPOINT
        return length(me.in1);
    }

    template < typename TIn1, typename TIn2, typename TIn3, typename TIn4, typename TIn5 >
    inline typename Size< Bundle5<TIn1, TIn2, TIn3, TIn4, TIn5> >::Type
    length(Bundle5<TIn1, TIn2, TIn3, TIn4, TIn5> const &me) {
SEQAN_CHECKPOINT
        return length(me.in1);
    }

//////////////////////////////////////////////////////////////////////////////


	template < typename TInput, typename TSpec >
    inline typename Size< Pipe<TInput, TSpec> >::Type
    countSequences(Pipe<TInput, TSpec> const &me) {
SEQAN_CHECKPOINT
        return countSequences(me.in);
    }

    template < typename TInput1, typename TInput2 >
    inline typename Size< Bundle2<TInput1, TInput2> >::Type
    countSequences(Bundle2<TInput1, TInput2> const &me) {
SEQAN_CHECKPOINT
        return countSequences(me.in1);
    }

    template < typename TInput1, typename TInput2, typename TInput3 >
    inline typename Size< Bundle3<TInput1, TInput2, TInput3> >::Type
    countSequences(Bundle3<TInput1, TInput2, TInput3> const &me) {
SEQAN_CHECKPOINT
        return countSequences(me.in1);
    }

    template < typename TIn1, typename TIn2, typename TIn3, typename TIn4, typename TIn5 >
    inline typename Size< Bundle5<TIn1, TIn2, TIn3, TIn4, TIn5> >::Type
    countSequences(Bundle5<TIn1, TIn2, TIn3, TIn4, TIn5> const &me) {
SEQAN_CHECKPOINT
        return countSequences(me.in1);
    }
/**
.Function.Pipe#front
..class:Class.Pipe
..cat:Pipelining
..summary:Gets the first element of the remaining stream.
..signature:front(object)
..param.object:A pop-passive pipeline module.
...type:Class.Pipe
..returns:The first element of the remaining input stream.
Return type is $Value<TObject>::Type$ for $object$ type $TObject$.
..remarks:@Function.Pipe#front@ or @Function.pop@ can only be called within a read process surrounded by @Function.beginRead@ and @Function.endRead@.
..see:Function.pop
..include:seqan/pipe.h
*/

    template < typename TInput, typename TSpec, typename TValue >
    inline typename Value< Pipe<TInput, TSpec> >::Type const &
    front(Pipe<TInput, TSpec> &me) {
SEQAN_CHECKPOINT
        return *me;
    }

/**
.Function.pop:
..class:Class.Pipe
..cat:Pipelining
..summary:Pops the first element of the remaining stream.
..signature:pop(object[, ref])
..param.object:A pop-passive pipeline module.
...type:Class.Pipe
..param.ref:Reference to the result. Result type is $Value<TObject>::Type$ for $object$ type $TObject$.
...remarks:Returns the first element of the remaining input stream.
..remarks:In contrast to @Function.Pipe#front@ this function also steps one element further.
..remarks:@Function.Pipe#front@ or @Function.pop@ can only be called within a read process surrounded by @Function.beginRead@ and @Function.endRead@.
..DISABLED.see:Function.top
..include:seqan/pipe.h
*/

    template < typename TInput, typename TSpec, typename TValue >
    inline void pop(Pipe<TInput, TSpec> &me, TValue &Ref_) {
SEQAN_CHECKPOINT
        Ref_ = *me;
        ++me;
    }

    template < typename TInput, typename TSpec >
    inline void pop(Pipe<TInput, TSpec> &me) {
SEQAN_CHECKPOINT
        ++me;
    }

///.Function.atEnd.param.iterator.type:Class.Pipe
///.Function.atEnd.class:Class.Pipe


    //////////////////////////////////////////////////////////////////////////////
    // pipe flow control

	struct ControlEof_;			// end of stream
	struct ControlEos_;			// end of sequence (for multiple sequences)
	struct ControlClear_;		// clear previous pool
	struct ControlBeginRead_;	// begin read process
	struct ControlEndRead_;		// end read process

	typedef Tag<ControlEof_>		ControlEof;
	typedef Tag<ControlEos_>		ControlEos;
	typedef Tag<ControlClear_>		ControlClear;
	typedef Tag<ControlBeginRead_>	ControlBeginRead;
	typedef Tag<ControlEndRead_>	ControlEndRead;

    template < typename TInput, typename TSpec, typename TCommand >
	inline bool control(Pipe<TInput, TSpec> &me, TCommand const &command) {
SEQAN_CHECKPOINT
        return control(me.in, command);
    }

    template < typename TInput, typename TSpec >
	inline bool eof(Pipe<TInput, TSpec> &me) {
SEQAN_CHECKPOINT
        return control(me, ControlEof());
    }

    template < typename TInput, typename TSpec >
	inline bool eos(Pipe<TInput, TSpec> &me) {
SEQAN_CHECKPOINT
        return control(me, ControlEos());
    }

    template < typename TInput, typename TSpec >
	inline bool clear(Pipe<TInput, TSpec> &me) {
SEQAN_CHECKPOINT
        return control(me, ControlClear());
    }

/**
.Function.beginRead
..class:Class.Pipe
..class:Class.Pool
..cat:Pipelining
..summary:Initiates a read process.
..signature:beginRead(object)
..param.object:A pop-passive pipeline module.
...type:Class.Pipe
...type:Class.Pool
..returns:A $bool$ which is $true$ on success.
..remarks:$beginRead$ rewinds the output stream, prepares $object$ for succeeding reads, and typically calls $beginRead$ of the input pipeline modules.
..remarks:A read process must be terminated with @Function.endRead@. Nested read processes are not allowed.
..see:Function.endRead
..include:seqan/pipe.h
*/

    template < typename TInput, typename TSpec >
	inline bool beginRead(Pipe<TInput, TSpec> &me) {
        return control(me, ControlBeginRead());
    }

/**
.Function.endRead
..class:Class.Pipe
..class:Class.Pool
..cat:Pipelining
..summary:Terminates a read process.
..signature:beginRead(object)
..param.object:A pop-passive pipeline module.
...type:Class.Pipe
...type:Class.Pool
..returns:A $bool$ which is $true$ on success.
..remarks:$endRead$ closes the output stream, frees resources possibly allocated by @Function.beginRead@, and typically calls $endRead$ of the input pipeline modules.
..see:Function.beginRead
..include:seqan/pipe.h
*/

    template < typename TInput, typename TSpec >
	inline bool endRead(Pipe<TInput, TSpec> &me) {
        return control(me, ControlEndRead());
    }


    //////////////////////////////////////////////////////////////////////////////
    // 2->1 pipe flow control
    template < typename TInput1, typename TInput2, typename TCommand >
    inline bool control(Bundle2<TInput1, TInput2> &me, TCommand const &command) {
        return	control(me.in1, command) &&
				control(me.in2, command);
    }

    //////////////////////////////////////////////////////////////////////////////
    // 3->1 pipe flow control
    template < typename TInput1, typename TInput2, typename TInput3, typename TCommand >
    inline bool control(Bundle3<TInput1, TInput2, TInput3> &me, TCommand const &command) {
        return	control(me.in1, command) &&
				control(me.in2, command) &&
				control(me.in3, command);
    }

    //////////////////////////////////////////////////////////////////////////////
    // 5->1 pipe flow control
    template < typename TIn1, typename TIn2, typename TIn3, typename TIn4, typename TIn5, typename TCommand >
    inline bool control(Bundle5<TIn1, TIn2, TIn3, TIn4, TIn5 > &me, TCommand const &command) {
        return	control(me.in1, command) &&
				control(me.in2, command) &&
				control(me.in3, command) &&
				control(me.in4, command) &&
				control(me.in5, command);
    }

///.Function.assign.param.source.type:Class.Pipe
///.Function.assign.class:Class.Pipe

	//////////////////////////////////////////////////////////////////////////////
    // pipe -> string

    // We cannot use the most-generic TObject as first argument, as operator << specialized somewhere else in
    // the first argument and here we need to specialize it in the second.
    template < typename TValue,
               typename TStringSpec,
               typename TInput,
               typename TSpec >
    inline void assign(String<TValue, TStringSpec> &dest, Pipe<TInput, TSpec> &src)
    {
        typedef typename Iterator<String<TValue, TStringSpec>, Standard >::Type TDestIter;
        resize(dest, length(src));
        beginRead(src);
        for (TDestIter _cur = begin(dest, Standard()), _end = end(dest, Standard()); _cur != _end; ++_cur, ++src)
            *_cur = *src;
        endRead(src);
    }

    template < typename TValue,
               typename TSegmentSpec,
               typename TInput,
               typename TSpec >
    inline void assign(Segment<TValue, TSegmentSpec> &dest, Pipe<TInput, TSpec> &src)
    {
        typedef typename Iterator<Segment<TValue, TSegmentSpec>, Standard >::Type TDestIter;
        resize(dest, length(src));
        beginRead(src);
        for (TDestIter _cur = begin(dest, Standard()), _end = end(dest, Standard()); _cur != _end; ++_cur, ++src)
            *_cur = *src;
        endRead(src);
    }

    template < typename TValue,
               typename TStringSpec,
               typename TInput,
               typename TSpec >
    inline void operator << (String<TValue, TStringSpec> &dest, Pipe<TInput, TSpec> &src)
    {
        assign(dest, src);
    }

    template < typename TValue,
               typename TSegmentSpec,
               typename TInput,
               typename TSpec >
    inline void operator << (Segment<TValue, TSegmentSpec> &dest, Pipe<TInput, TSpec> &src)
    {
        assign(dest, src);
    }

	//////////////////////////////////////////////////////////////////////////////
    // pipe -> out_stream
    template <typename TInput, typename TSpec>
	std::ostream& operator << (std::ostream &out, Pipe<TInput, TSpec> &src)
    {
        beginRead(src);
        while (!eof(src)) {
			out << *src << ::std::endl;
            ++src;
        }
        endRead(src);
		return out;
	}


    template < typename TObject, typename TSpec >
    struct BufferHandler;

    template <typename TObject, typename TSpec>
    struct Value<BufferHandler<TObject, TSpec> >
    {
        typedef Buffer<typename Value<TObject>::Type> Type;
    };



    template <typename TObject, typename TSpec>
    struct Handler;

    template <typename TObject, typename TSpec>
    struct Value<Handler<TObject, TSpec> >:
        public Value<TObject> {};

    // buffer-based read/write handler metafunctions
    template < typename TInput >
    struct BufReadHandler;

    template < typename TOutput >
    struct BufWriteHandler;



	//////////////////////////////////////////////////////////////////////////////
	// generic adapter for buffered readers/writers

	template <typename TValue, typename TSpec, typename TSize >
    inline void resize(Buffer<TValue, TSpec> &me, TSize size);

	struct AdapterSpec;

	template <typename TBufferHandler>
    struct Handler<TBufferHandler, AdapterSpec>
    {
        typedef typename Value<TBufferHandler>::Type        TBuffer;
        typedef typename Value<TBuffer>::Type               TValue;
		typedef typename Iterator<TBuffer, Standard>::Type  TIterator;

        TBufferHandler  handler;
        TBuffer			buffer;
        TIterator       cur;

        template < typename TObject >
        Handler(TObject &_object):
            handler(_object) {}

        inline bool begin()
        {
            buffer = handler.first();
            cur = seqan::begin(buffer, Standard());
            return true;
        }

        inline TValue const & front() const
        {
            return *cur;
        }

        inline void pop()
        {
			if (++cur == seqan::end(buffer, Standard()))
            {
                buffer = handler.next();
				cur = seqan::begin(buffer, Standard());
			}
        }

        inline void pop(TValue &Ref_)
        {
            Ref_ = *cur;
            pop();
        }

        inline void push(TValue const & Val_)
        {
            if (cur == seqan::end(buffer, Standard()))
            {
                buffer = handler.next();
                cur = seqan::begin(buffer, Standard());
            }
            *cur = Val_;
            ++cur;
        }

        inline bool eof() const
        {
            return length(buffer) == 0;
        }

        inline void end()
        {
            handler.end();
            resize(buffer, 0);
        }

        inline void process()
        {
            handler.process();
        }
    };

    template <typename TBufferHandler>
    struct Value<Handler<TBufferHandler, AdapterSpec> >:
        public Value<typename Value<TBufferHandler>::Type> {};

    // character-based read/write handler metafunctions
    template < typename TInput >
    struct ReadHandler
    {
        typedef Handler< typename BufReadHandler<TInput> ::Type, AdapterSpec > Type;
    };

    template < typename TOutput >
    struct WriteHandler
    {
        typedef Handler< typename BufWriteHandler<TOutput> ::Type, AdapterSpec > Type;
    };


	//////////////////////////////////////////////////////////////////////////////
	// pair incrementer
	//
	// used by pipes processing multiples sequences 
	// for generating pairs (seqNo, seqOffs)

	template <typename TPair, typename TLimits>
	struct PairIncrementer_
    {
        typedef typename Iterator<TLimits const, Standard>::Type            TIter;
        typedef typename RemoveConst_<typename Value<TLimits>::Type>::Type  TSize;
        typedef typename Value<TPair, 2>::Type                              TOffset;

        TIter   it, itEnd;
		TSize   old;
        TOffset localEnd;

		TPair pos;
		inline operator TPair () const
        {
			return pos;
		}

		inline TPair const & operator++ ()
        {
			TOffset i2 = getValueI2(pos) + 1;
			if (i2 >= localEnd)
            {
				i2 = 0;
				localEnd = 0;
				while (!localEnd && (it != itEnd))
				{
                    typename Value<TPair,1>::Type i1 = getValueI1(pos);
                    assignValueI1(pos, i1 + 1);
					localEnd = (*it - old);
                    // test for overflows
                    SEQAN_ASSERT_LT_MSG(i1, getValueI1(pos), "Overflow detected. Use a bigger type for the *first* value in the SAValue pair!");
                    SEQAN_ASSERT_EQ_MSG((TSize)localEnd, (*it - old), "Overflow detected. Use a bigger type for the *second* value in the SAValue pair!");
					old = *it;
					++it;
				}
				if (!localEnd && it == itEnd)
					assignValueI1(pos, getValueI1(pos) + 1);	// set pos behind the last sequence
			}
			assignValueI2(pos, i2);
			return pos;
		}
	};

	template <typename TPair, typename TLimits, typename TLimits2>
	void setHost(PairIncrementer_<TPair, TLimits> &me, TLimits2 const &limits)
    {
		me.it = begin(limits, Standard());
		me.itEnd = end(limits, Standard());
		me.old = 0;
		me.localEnd = 0;
        me.pos = TPair(0, 0);
		if (length(limits) > 1)
        {
			++me.it;
			++me;
			assignValueI1(me.pos, getValueI1(me.pos) - 1);
		}
	}

	template <typename TPair, typename TLimits, typename TLimits2>
	void setHost(PairIncrementer_<TPair, TLimits> &me, TLimits2 &limits)
    {
        setHost(me, const_cast<TLimits2 const &>(limits));
    }
//____________________________________________________________________________

	template <typename TPair, typename TLimits>
	TPair const & value(PairIncrementer_<TPair, TLimits> const &me) {
		return me.pos;
	}

	template <typename TPair, typename TLimits>
	TPair & value(PairIncrementer_<TPair, TLimits> &me) {
		return me.pos;
	}


	//////////////////////////////////////////////////////////////////////////////
	// pair decrementer
	//
	// used by pipes processing multiples sequences 
	// for generating pairs (seqNo, seqOffs)

	template <typename TPair, typename TLimits, unsigned m = 0>
	struct PairDecrementer_
    {
        typedef typename Iterator<TLimits const, Standard>::Type            TIter;
        typedef typename RemoveConst_<typename Value<TLimits>::Type>::Type  TSize;
        
        TIter       it, itEnd;
		TSize       old;
		TPair		pos;
		unsigned	residue;

		PairDecrementer_() {}
		PairDecrementer_(TLimits const &_limits) { setHost(*this, _limits); }

		inline operator TPair () const {
			return pos;
		}

		inline TPair const & operator-- ()
        {
			typename Value<TPair,2>::Type i2 = getValueI2(pos);
			if (i2 > 1)
            {
				--i2;
				if (residue == 0) residue = m;
				--residue;
			}
			else
			{
				i2 = 0;
				while (!i2 && (it != itEnd))
				{
                    typename Value<TPair,1>::Type i1 = getValueI1(pos);
					assignValueI1(pos, i1 + 1);
					i2 = (*it - old);
                    // test for overflows
                    SEQAN_ASSERT_LT_MSG(i1, getValueI1(pos), "Overflow detected. Use a bigger type for the *first* value in the SAValue pair!");
                    SEQAN_ASSERT_EQ_MSG((TSize)i2, *it - old, "Overflow detected. Use a bigger type for the *second* value in the SAValue pair!");
					old = *it;
					++it;
				} 
				residue = i2 % m;
			}
			assignValueI2(pos, i2);
			return pos;
		}
	};

	template <typename TPair, typename TLimits, unsigned m, typename TLimits2>
	void setHost(PairDecrementer_<TPair, TLimits, m> &me, TLimits2 const &limits)
    {
		me.it = begin(limits);
		me.itEnd = end(limits);
		me.old = 0;
        me.pos = TPair(0, 0);
		if (length(limits) > 1)
        {
			++me.it;
			--me;
			assignValueI1(me.pos, getValueI1(me.pos) - 1);
		} else
			me.residue = 0;
	}

	template <typename TPair, typename TLimits, unsigned m, typename TLimits2>
	void setHost(PairDecrementer_<TPair, TLimits, m> &me, TLimits2 &limits)
    {
        setHost(me, const_cast<TLimits const &>(limits));
    }
//____________________________________________________________________________

	template <typename TPair, typename TLimits>
	struct PairDecrementer_<TPair, TLimits, 0>
    {
        typedef typename Iterator<TLimits const, Standard>::Type            TIter;
        typedef typename RemoveConst_<typename Value<TLimits>::Type>::Type  TSize;

        TIter       it, itEnd;
		TSize       old;
		TPair		pos;

		PairDecrementer_() {}
		PairDecrementer_(TLimits const &_limits) { setHost(*this, _limits); }

		inline operator TPair () const {
			return pos;
		}

		inline TPair const & operator-- ()
        {
			typename Value<TPair,2>::Type i2 = getValueI2(pos);
			if (i2 > 1)
				--i2;
			else
			{
				i2 = 0;
				while (!i2 && (it != itEnd))
				{
                    typename Value<TPair,1>::Type i1 = getValueI1(pos);
					assignValueI1(pos, i1 + 1);
					i2 = (*it - old);
                    // test for overflows
                    SEQAN_ASSERT_LT_MSG(i1, getValueI1(pos), "Overflow detected. Use a bigger type for the *first* value in the SAValue pair!");
                    SEQAN_ASSERT_EQ_MSG((TSize)i2, *it - old, "Overflow detected. Use a bigger type for the *second* value in the SAValue pair!");
					old = *it;
					++it;
				} 
			}
			assignValueI2(pos, i2);
			return pos;
		}
	};

	template <typename TPair, typename TLimits, typename TLimits2>
	void setHost(PairDecrementer_<TPair, TLimits, 0> &me, TLimits2 const &limits) {
		me.it = begin(limits);
		me.itEnd = end(limits);
		me.old = 0;
		assignValueI1(me.pos, 0);
		assignValueI2(me.pos, 0);
		if (length(limits) > 1) {
			++me.it;
			--me;
			assignValueI1(me.pos, getValueI1(me.pos) - 1);
		}
	}

	template <typename TPair, typename TLimits, typename TLimits2>
	void setHost(PairDecrementer_<TPair, TLimits, 0> &me, TLimits2 &limits)
    {
        setHost(me, const_cast<TLimits const &>(limits));
    }
//____________________________________________________________________________

	template <typename TPair, typename TLimits, unsigned m>
	TPair const & value(PairDecrementer_<TPair, TLimits, m> const &me) {
		return me.pos;
	}

	template <typename TPair, typename TLimits, unsigned m>
	TPair & value(PairDecrementer_<TPair, TLimits, m> &me) {
		return me.pos;
	}

}

#endif
