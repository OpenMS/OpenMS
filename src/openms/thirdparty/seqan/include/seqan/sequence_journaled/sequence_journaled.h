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
// Journaled String implementation.
// ==========================================================================

// TODO(holtgrew): Journaled strings will probably not work for non-POD alphabets!

#ifndef SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_
#define SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Specialization Journaled String
// ----------------------------------------------------------------------------

/**
.Spec.Journaled String
..cat:Sequences
..general:Class.String
..summary:Journaled versions of arbitrary underlying string.
..signature:String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >
..include:seqan/sequence_journaled.h
 */

template <typename THostSpec, typename TJournalSpec = SortedArray, typename TBufferSpec = Alloc<void> >
struct Journaled {};

template <typename TValue_, typename THostSpec_, typename TJournalSpec_, typename TBufferSpec_>
class String<TValue_, Journaled<THostSpec_, TJournalSpec_, TBufferSpec_> >
{
public:
    typedef String<TValue_, Journaled<THostSpec_, TJournalSpec_, TBufferSpec_> > TThis_;
    
    typedef TValue_ TValue;
    typedef THostSpec_ THostSpec;
    typedef TJournalSpec_ TJournalSpec;
    typedef TBufferSpec_ TBufferSpec;

    typedef String<TValue, THostSpec> THost;
    typedef typename Size<THost>::Type TSize;
    typedef typename Position<THost>::Type TPosition;
    typedef String<TValue, TBufferSpec> TInsertionBuffer;
    typedef JournalEntry<TSize, TPosition> TJournalEntry;
    typedef JournalEntries<TJournalEntry, TJournalSpec> TJournalEntries;

    // The underlying host string.
    Holder<THost> _holder;
    // A buffer for inserted strings.
    TInsertionBuffer _insertionBuffer;
    // The journal is a sorted set of TJournalEntry objects, the exact types
    // depends on TJournalSpec.  Note that the entries resemble a partial
    // sum datastructure.
    TJournalEntries _journalEntries;
    // The journaled string's size.
    TSize _length;

    String() {
        SEQAN_CHECKPOINT;
    }

    // Note: Defining both, constructors from same type and other for clarity.
    
    String(THost & host)
    {
        SEQAN_CHECKPOINT;
        setHost(*this, host);
    }

    String(String const & other)
    {
        SEQAN_CHECKPOINT;
        assign(_holder, other._holder);
        set(*this, other);
    }

    template <typename TString>
	String(TString const & other)
    {
        SEQAN_CHECKPOINT;
        assign(*this, other);
    }

    inline
    String &
    operator=(String const & other)
    {
        SEQAN_CHECKPOINT;
        //TODO (rmaerker): should copy the holder, too. What is intended with the assignment - copy or set?
        //TODO (rmaerker): missing self-assignment test. What happens with ptr of Holder in self-assignment
        // assign(_holder, other._holder);
        assign(*this, other);
        return *this;
    }

    template <typename TString>
    inline
    String &
    operator=(TString const & other)
    {
        SEQAN_CHECKPOINT;
        assign(*this, other);
        return *this;
    }

    inline
    typename Reference<String>::Type
    operator[](TPosition pos)
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }

    inline
    typename Reference<String const>::Type
    operator[](TPosition pos) const
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

/**
.Metafunction.Host
..param.T.type:Spec.Journaled String
..class:Spec.Journaled String
*/

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef String<TValue, THostSpec> Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
    typedef String<TValue, THostSpec> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction InsertionBuffer
// ----------------------------------------------------------------------------

/**
.Metafunction.InsertionBuffer
..cat:Sequences
..summary:Return type of insertion buffer string for a journaled string.
..class:Spec.Journaled String
..signature:InsertionBuffer<T>::Type
..param.T:The journaled string to get the insertion buffer type for.
...type:Spec.Journaled String
..include:sequan/sequence_journal.h
*/
template <typename T>
struct InsertionBuffer;

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct InsertionBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef String<TValue, TBufferSpec> Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct InsertionBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
    typedef String<TValue, TBufferSpec> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

///.Metafunction.Size.param.T.type:Spec.Journaled String
///.Metafunction.Size.class:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef typename String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >::TSize Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
    : Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > > {};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

///.Metafunction.Position.param.T.type:Spec.Journaled String
///.Metafunction.Position.class:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
  typedef typename String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >::TPosition Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
    : Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > > {};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

///.Metafunction.Reference.param.T.type:Spec.Journaled String
///.Metafunction.Reference.class:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Reference<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > TJournaledString;
    typedef typename Iterator<TJournaledString>::Type TIterator_;
    typedef Proxy<IteratorProxy<TIterator_> > TProxy_;
    typedef TProxy_ Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Reference<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const TJournaledString;
    typedef typename Iterator<TJournaledString>::Type TIterator_;
    typedef Proxy<IteratorProxy<TIterator_> > TProxy_;
    typedef TProxy_ Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

///.Metafunction.GetValue.param.T.type:Spec.Journaled String
///.Metafunction.GetValue.class:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct GetValue<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
  typedef TValue Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct GetValue<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
  typedef TValue Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// TODO(holtgrew): Does Value have to be overwritten? Is not for packed strings!
///.Metafunction.Value.param.T.type:Spec.Journaled String
///.Metafunction.Value.class:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Value<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
  typedef TValue Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Value<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
  typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction JournalType
// ----------------------------------------------------------------------------

/**
.Metafunction.JournalType
..cat:Sequences
..summary:Metafunction for returning the type of the journal of a Journaled String.
..signature:JournalType<T>::Type
..class:Spec.Journaled String
..param.T:The journaled string to get the journal type for.
...type:Spec.Journaled String
..include:seqan/string_journaled.h
 */
template <typename T>
struct JournalType;

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct JournalType<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef typename Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type TSize_;
    typedef typename Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type TPosition_;
    typedef JournalEntry<TSize_, TPosition_> TJournalEntry_;

    typedef JournalEntries<TJournalEntry_, TJournalSpec> Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct JournalType<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
    typedef typename JournalType<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator<<
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
TStream &
operator<<(TStream & stream, String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & s)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > TString;
    typedef typename TString::TJournalEntries TJournalEntries;
    typedef typename Iterator<TJournalEntries const, Standard>::Type TIterator;
    
    for (TIterator it = begin(s._journalEntries), itend = end(s._journalEntries); it != itend; ++it) {
        if (value(it).segmentSource == SOURCE_ORIGINAL) {
            stream << infix(value(s._holder), value(it).physicalPosition, value(it).physicalPosition + value(it).length);
        } else {
            SEQAN_ASSERT_EQ(value(it).segmentSource, SOURCE_PATCH);
            stream << infix(s._insertionBuffer, value(it).physicalPosition, value(it).physicalPosition + value(it).length);
        }
    }
    return stream;
}

// ----------------------------------------------------------------------------
// Function assign
// ----------------------------------------------------------------------------

// Assignment always resizes the insertion buffer and copies the contents of
// source into the insertion buffer.  Even when assigning a journaled string
// to a journaled string!
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TSource>
inline
void
assign(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    clear(target);
    replace(target, 0, length(target), source);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TSource>
inline
void
assign(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, static_cast<TSource const &>(source));
}

// ----------------------------------------------------------------------------
// Function set
// ----------------------------------------------------------------------------

// Setting a journaled string to another journaled string copies over the deep
// structure.  For all other cases, it is assign().
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
void
set(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & target,
    String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & source)
{
    SEQAN_CHECKPOINT;
    target._holder = source._holder;
    target._insertionBuffer = source._insertionBuffer;
    target._journalEntries = source._journalEntries;
    target._length = source._length;
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
void
set(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & target,
    String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & source)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > TJournaledString;
    set(target, static_cast<TJournaledString const &>(source));
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TSource>
inline
void
set(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & target,
    TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

// ----------------------------------------------------------------------------
// Function setHost
// ----------------------------------------------------------------------------

/**
.Function.setHost:
..class:Spec.Journaled String
..param.object.type:Spec.Journaled String
..param.host.type:Class.String
..include:seqan/sequence_journaled.h
*/
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TSequence2>
inline
void
setHost(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString, TSequence2 & str)
{
    SEQAN_CHECKPOINT;
    setValue(journaledString._holder, str);
    journaledString._length = length(str);
    reinit(journaledString._journalEntries, length(str));
}

// ----------------------------------------------------------------------------
// Function host
// ----------------------------------------------------------------------------

/**
.Function.host:
..class:Spec.Journaled String
..param.object.type:Spec.Journaled String
..include:seqan/sequence_journaled.h
*/
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type &
host(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString)
{
    SEQAN_CHECKPOINT;
    return value(journaledString._holder);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type const &
host(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString)
{
    SEQAN_CHECKPOINT;
    return value(journaledString._holder);
}

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

/**
.Function.clear:
..class:Spec.Journaled String
..param.object.type:Spec.Journaled String
..include:seqan/sequence_journaled.h
 */
// TODO(holtgrew): Behaviour is to clear the journal, not the string!
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline void
clear(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString)
{
    SEQAN_CHECKPOINT;
    reinit(journaledString._journalEntries, length(host(journaledString)));
    clear(journaledString._insertionBuffer);
    _setLength(journaledString, length(host(journaledString)));
}

// ----------------------------------------------------------------------------
// Function flatten
// ----------------------------------------------------------------------------

/**
.Function.flatten:
..class:Spec.Journaled String  
..cat:Sequences
..summary:Apply the journal to the underlying string, destructively on the underlying string.
..signature:flatten(journaledString)
..param.journaledString:The journaled string to flatten.
...type:Spec.Journaled String
..include:seqan/sequence_journaled.h
 */
// TODO(holtgrew): What about non-destructive version that creates a new copy and sets holder to it?
// TODO(rmaerker): Only supported by SortedArray, since there is a bug in the UnbalancedTree iterator. Is not fixed since UT is deprecated.
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline void
flatten(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > TJournalString;
    typedef typename Host<TJournalString>::Type THost;
    typedef typename Position<TJournalString>::Type TPosition;
    typedef typename JournalType<TJournalString>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TEntriesIterator;

    TEntriesIterator it = end(journaledString._journalEntries);
    TEntriesIterator itBegin = begin(journaledString._journalEntries);

    TPosition lastRefPos = length(host(journaledString));
    while (it != itBegin)
    {
        --it;
        if (value(it).segmentSource == SOURCE_ORIGINAL)
        {
            // Check for a deletion and if there is one then erase this part from reference.
            if ((value(it).physicalPosition + value(it).length) < lastRefPos)
                erase(host(journaledString), value(it).physicalPosition + value(it).length, lastRefPos);
            lastRefPos = value(it).physicalPosition;
        }
        else if (value(it).segmentSource == SOURCE_PATCH)
        {
            insert(host(journaledString),
                   lastRefPos,
                   THost(infix(journaledString._insertionBuffer, value(it).physicalPosition,
                               value(it).physicalPosition + value(it).length)));
        }
    }
    // Handle deletion at beginning of reference.
    SEQAN_ASSERT_GEQ(lastRefPos, (TPosition) 0);
    if (0 != lastRefPos)
        erase(host(journaledString), (TPosition) 0, lastRefPos);
    // After transmitting the journaled differences to the reference clear the insertion buffer and reset everything.
    clear(journaledString);
}


// ----------------------------------------------------------------------------
// Function erase
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TBeginPos, typename TEndPos>
inline void
erase(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
      TBeginPos pos,
      TEndPos posEnd)
{
    SEQAN_CHECKPOINT;
	SEQAN_ASSERT_GEQ(static_cast<TBeginPos>(journaledString._length), pos);
	SEQAN_ASSERT_GEQ(static_cast<TEndPos>(journaledString._length), posEnd);
    SEQAN_ASSERT_GEQ(static_cast<TBeginPos>(journaledString._length), static_cast<TBeginPos>(posEnd - pos));
    journaledString._length -= posEnd - pos;
    recordErase(journaledString._journalEntries, pos, posEnd);
    if (length(journaledString._journalEntries) == 0)
        clear(journaledString._insertionBuffer);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline void
erase(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(journaledString._length, 1u);
    erase(journaledString, pos, pos + 1);
}

// ----------------------------------------------------------------------------
// Function insert
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TString, typename TPos>
inline void
insert(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
       TPos pos,
       TString const & seq)
{
    SEQAN_CHECKPOINT;
    journaledString._length += length(seq);
    TPos beginPos = length(journaledString._insertionBuffer);
    append(journaledString._insertionBuffer, seq);
    recordInsertion(journaledString._journalEntries, pos, beginPos, length(seq));
}

// ----------------------------------------------------------------------------
// Function insertValue
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos, typename TValue2>
inline void
insertValue(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
            TPos pos,
            TValue2 const & value)
{
    SEQAN_CHECKPOINT;
    journaledString._length += 1;
    TPos beginPos = length(journaledString._insertionBuffer);
    appendValue(journaledString._insertionBuffer, value);
    recordInsertion(journaledString._journalEntries, pos, beginPos, 1u);
}

// ----------------------------------------------------------------------------
// Function assignInfix
// ----------------------------------------------------------------------------

// TODO(holtgrew): More of an undescore function.

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TBeginPos, typename TEndPos, typename TSequence2>
inline void
assignInfix(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
            TBeginPos beginPos,
            TEndPos endPos,
            TSequence2 const & valueString)
{
    SEQAN_CHECKPOINT;
    erase(journaledString, beginPos, endPos);
    insert(journaledString, beginPos, valueString);
}

// ----------------------------------------------------------------------------
// Function assignValue
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos, typename TValue2>
inline void
assignValue(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
            TPos pos,
            TValue2 const & value)
{
    SEQAN_CHECKPOINT;
    erase(journaledString, pos);
    insertValue(journaledString, pos, value);
}


// TODO(holtgrew): Batch-Assignment of values through segments?

// TODO(holtgrew): begin
// TODO(holtgrew): empty
// TODO(holtgrew): end
// TODO(holtgrew): flatten
// TODO(holtgrew): fill
// TODO(holtgrew): getValue

// TODO(holtgrew): Unused, remove?
/*
template <typename TSequence, typename TJournalSpec>
inline
typename Value<TSequence>::Type const &
front(String...<TSequence, TJournalSpec> const & journaledString)
{
    SEQAN_XXXCHECKPOINT;
    typedef SequenceJournal<TSequence, TJournalSpec> TString;
    typedef typename TString::TNode TNode;
    TNode frontNode = front(journaledString._journalEntries);
    if (frontNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journaledString._holder), frontNode->virtualPosition + frontNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(frontNode->segmentSource, SOURCE_PATCH);
        return getValue(journaledString._insertionBuffer, frontNode->virtualPosition + frontNode->length - 1);
    }
}

// front/back clash with general sequence definitions.
template <typename TSequence, typename TJournalSpec>
inline
TValue const &
back(SequenceJournal<TSequence, TJournalSpec> const & journaledString)
{
    SEQAN_XXXCHECKPOINT;
    typedef SequenceJournal<TSequence, TJournalSpec> TString;
    typedef typename TString::TNode TNode;
    TNode backNode = back(journaledString._journalEntries);
    if (backNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journaledString._holder), backNode->virtualPosition + backNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(backNode->segmentSource, SOURCE_PATCH);
        return getValue(journaledString._insertionBuffer, backNode->virtualPosition + backNode->length - 1);
    }
}
*/

// ----------------------------------------------------------------------------
// Function length
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type
length(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString)
{
    SEQAN_CHECKPOINT;
    return journaledString._length;
}

// ----------------------------------------------------------------------------
// Function toCString
// ----------------------------------------------------------------------------

// TODO(holtgrew): toCString

// ----------------------------------------------------------------------------
// Function value
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline
typename Reference<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type
value(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & me, 
      TPos pos)
{
    SEQAN_CHECKPOINT;
    
    return *iter(me, pos, Standard());
} 

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline
typename Reference<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>::Type
value(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & me, 
      TPos pos)
{
    SEQAN_CHECKPOINT;
    
    return *iter(me, pos, Standard());
} 


// TODO(holtgrew): Maybe better use template parameter TPos?
// TOOD(holtgrew): operator<
// TOOD(holtgrew): operator>
// TOOD(holtgrew): operator<=
// TOOD(holtgrew): operator>=
// TOOD(holtgrew): operator==
// TOOD(holtgrew): operator!=

// ----------------------------------------------------------------------------
// Function getValue
// ----------------------------------------------------------------------------

// TODO(holtgrew): Maybe better use template parameter TPos?
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename GetValue<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>::Type
getValue(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString,
         typename Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type pos)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const TJournaledString;
    typedef typename TJournaledString::TJournalEntry TJournalEntry;
    typedef typename Position<TJournaledString>::Type TPos;

    TJournalEntry entry = findJournalEntry(journaledString._journalEntries, pos);
    TPos relativePos = pos - entry.virtualPosition;

    if (entry.segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journaledString._holder), entry.physicalPosition + relativePos);
    } else {
        return getValue(journaledString._insertionBuffer, entry.physicalPosition + relativePos);
    }
}

// --------------------------------------------------------------------------
// Function virtualToHostPosition()
// --------------------------------------------------------------------------

/**
.Function.Journaled String#virtualToHostPosition
..class:Spec.Journaled String
..summary:Translates virtual (view) position to position in host.
..signature:virtualToHostPosition(journaledString, pos)
..param.journaledString:The journaled string to translate for.
...type:Spec.Journaled String
..param.pos:The position to translate.
..returns:Corresponding host position.
...type:Metafunction.Position
*/

// Note that if pos is in a gap, we return the position of the entry
// after the gap in the host.
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline
typename Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type
virtualToHostPosition(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString,
                      TPos pos)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): With a better journal entries datastructure, we could solve the main problem here. At the moment, we delegate completely.
    return virtualToHostPosition(journaledString._journalEntries, pos);
}

//---------------------------------------------------
// function hostToVirtualPosition()
//---------------------------------------------------

/**
.Function.Journaled String#hostToVirtualPosition
..class:Spec.Journaled String
..summary:Translates position in host to virtual (view) position.
..signature:hostToVirtualPosition(journaledString, pos)
..param.journaledString:The journaled string to translate for.
...type:Spec.Journaled String
..param.pos:The position to translate.
..returns:Corresponding virtual (view) position.
...type:Metafunction.Position
*/

template<typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline
typename Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type
hostToVirtualPosition(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journalString,
                      TPos const & pos)
{
    SEQAN_ASSERT_LEQ(pos, (TPos) length(host(journalString)));
    return hostToVirtualPosition(journalString._journalEntries, pos);
}

// --------------------------------------------------------------------------
// Function isGapInHost()
// --------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline
bool
isGapInHost(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString,
            TPos pos)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): With a better journal entries datastructure, we could solve the main problem here. At the moment, we delegate completely.
    return isGapInHost(journaledString._journalEntries, pos);
}

// --------------------------------------------------------------------------
// Function _setLength()
// --------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
void
_setLength(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
           size_t newLength)
{
    SEQAN_CHECKPOINT;
    journaledString._length = newLength;
}

// --------------------------------------------------------------------------
// Function replace()
// --------------------------------------------------------------------------

template <typename TTargetValue, typename TTargetHostSpec, typename TTargetJournalSpec, typename TTargetBufferSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(String<TTargetValue, Journaled<TTargetHostSpec, TTargetJournalSpec, TTargetBufferSpec> > & target,
        TPositionBegin posBegin,
        TPositionEnd posEnd,
        TSource const & source,
        Tag<TExpand> /*tag*/)
{
    SEQAN_CHECKPOINT;
    assignInfix(target, posBegin, posEnd, source);
}

template <typename TTargetValue, typename TTargetHostSpec, typename TTargetJournalSpec, typename TTargetBufferSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(String<TTargetValue, Journaled<TTargetHostSpec, TTargetJournalSpec, TTargetBufferSpec> > & target,
        TPositionBegin posBegin,
        TPositionEnd posEnd,
        TSource const & source,
        typename Size<String<TTargetValue, Journaled<TTargetHostSpec, TTargetJournalSpec, TTargetBufferSpec> > >::Type limit,
        Tag<TExpand> tag)
{
    // Possibly only shorten target if limit is too small.
    if (posBegin > static_cast<TPositionBegin>(limit)) {
        assignInfix(target, limit, length(target), infix(source, 0, 0));
        return;
    }

    // Replace range in target and afterwards, limit length of target.
    replace(target, posBegin, posEnd, source, tag);
    assignInfix(target, limit, length(target), infix(source, 0, 0));
}

// This variant is a workaround for the "const array"-bug of VC++.
template <typename TTargetValue, typename TTargetHostSpec, typename TTargetJournalSpec, typename TTargetBufferSpec, typename TPositionBegin, typename TPositionEnd, typename TSourceValue, typename TExpand>
inline void
replace(String<TTargetValue, Journaled<TTargetHostSpec, TTargetJournalSpec, TTargetBufferSpec> > & target,
        TPositionBegin posBegin,
        TPositionEnd posEnd,
        TSourceValue const * source,
        Tag<TExpand> /*tag*/)
{
    SEQAN_CHECKPOINT;
    assignInfix(target, posBegin, posEnd, source);
}

// This variant is a workaround for the "const array"-bug of VC++.
template <typename TTargetValue, typename TTargetHostSpec, typename TTargetJournalSpec, typename TTargetBufferSpec, typename TPositionBegin, typename TPositionEnd, typename TSourceValue, typename TExpand>
inline void
replace(String<TTargetValue, Journaled<TTargetHostSpec, TTargetJournalSpec, TTargetBufferSpec> > & target,
        TPositionBegin posBegin,
        TPositionEnd posEnd,
        TSourceValue const * source,
        typename Size<String<TTargetValue, Journaled<TTargetHostSpec, TTargetJournalSpec, TTargetBufferSpec> > >::Type limit,
        Tag<TExpand> tag)
{
    // Possibly only shorten target if limit is too small.
    if (posBegin > static_cast<TPositionBegin>(limit)) {
        assignInfix(target, limit, length(target), infix(source, 0, 0));
        return;
    }

    // Replace range in target and afterwards, limit length of target.
    replace(target, posBegin, posEnd, source, tag);
    assignInfix(target, limit, length(target), infix(source, 0, 0));
}

// --------------------------------------------------------------------------
// Function getObjectId()
// --------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
const void *
getObjectId(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString)
{
    SEQAN_CHECKPOINT;
    return getObjectId(value(journaledString._holder));
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
const void *
getObjectId(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString)
{
    SEQAN_CHECKPOINT;
    return getObjectId(value(journaledString._holder));
}

// --------------------------------------------------------------------------
// Function isFlat()
// --------------------------------------------------------------------------

/**
.Function.Journaled String#isFlat
..class:Spec.Journaled String
..cat:Sequences
..summary:Returns whether journaled string has modifications.
..signature:isFlat(journaledString)
..param.journaledString:The journaled string to query.
...type:Spec.Journaled String
..returns:$bool$ indicating whether the string has been modfieid.
*/

// TODO(holtgrew): Is the non-const version necessary?
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
inline bool
isFlat(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString)
{
	SEQAN_CHECKPOINT;
	typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
	typedef typename JournalType<TJournalString>::Type TJournalEntries;
	typedef typename Iterator<TJournalEntries const, Standard>::Type TIteraror;

	TIteraror it = begin(journaledString._journalEntries);
	if ((*it).segmentSource == SOURCE_ORIGINAL)
	{
		if (((*it).physicalPosition == (*it).virtualPosition) && ((*it).length == length(host(journaledString))))
		{
			return true;
		}
	}
	return false;
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
inline bool
isFlat(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & journaledString)
{
	SEQAN_CHECKPOINT;
	typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
	typedef typename JournalType<TJournalString>::Type TJournalEntries;
	typedef typename Iterator<TJournalEntries const, Standard>::Type TIteraror;

	TIteraror it = begin(journaledString._journalEntries);
	if ((*it).segmentSource == SOURCE_ORIGINAL)
	{
		if (((*it).physicalPosition == (*it).virtualPosition) && ((*it).length == length(host(journaledString))))
		{
			return true;
		}
	}
	return false;
}

// template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
// void
// _deallocateStorage(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & string,
//				   size_t const & newCapacity)
// {
//	std::cout << "Now he is here but does not actually delete the contents." << std::endl;
// }

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_
