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
// Code for the Journaled string iterator.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
#define SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

template <typename TJournaledStringSpec>
struct JournaledStringIterSpec;

template <typename TJournaledString, typename TJournalSpec>
class Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
{
public:
    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;
    typedef typename TJournaledString::TValue TValue;
    typedef typename JournalType<TJournaledString>::Type TJournalEntries;
    // We need a rooted iterator for iterating the journal tree since we need atEnd().
    typedef typename Iterator<TJournalEntries, Rooted>::Type TJournalEntriesIterator;
    typedef typename Host<TJournaledString>::Type THost;
    typedef typename Iterator<THost, Standard>::Type THostIterator;
    typedef typename InsertionBuffer<TJournaledString>::Type TInsertionBuffer;
    typedef typename Iterator<TInsertionBuffer, Standard>::Type TInsertionBufferIterator;

    // The journal string we iterate over.
    TJournaledString * _journalStringPtr;
    // Iterator over the segments in the journal tree.
    TJournalEntriesIterator _journalEntriesIterator;
    // Begin and end iterator in the host string of the journal string.
    THostIterator _hostSegmentBegin;
    THostIterator _hostSegmentEnd;
    // Current iterator in the host segment.
    THostIterator _currentHostIt;
    // Begin and end iterator in the insertion buffer of the journal string.
    TInsertionBufferIterator _insertionBufferSegmentBegin;
    TInsertionBufferIterator _insertionBufferSegmentEnd;
    // Current iterator in the insertion buffer.
    TInsertionBufferIterator _currentInsertionBufferIt;

    Iter() :
        _journalStringPtr(), _journalEntriesIterator(), _hostSegmentBegin(), _hostSegmentEnd(),
        _currentHostIt(), _insertionBufferSegmentBegin(), _insertionBufferSegmentEnd(),
        _currentInsertionBufferIt()
    {}

    Iter(TIterator const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    {
        SEQAN_CHECKPOINT;
    }

    Iter(typename IterComplementConst<TIterator>::Type const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    {
        SEQAN_CHECKPOINT;
    }

    // TODO(holtgrew): Commented out because of weird ambiguities with other constructor.
    // explicit
    // Iter(TJournaledString & journalString)
    // {
    //     SEQAN_CHECKPOINT;
    //     _initJournaledStringIterator(*this, journalString);
    // }
};

// ============================================================================
// Metafunctions
// ============================================================================

// For String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >

///.Metafunction.Iterator.param.T:Spec.Journal String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>
{
    typedef Iter<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, JournaledStringIterSpec<TJournalSpec> > Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>
{
    typedef Iter<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, JournaledStringIterSpec<TJournalSpec> > Type;
};

// For Iter<TJournaledString, TJournaledStringIterSpec>
template <typename TJournaledString, typename TJournaledStringIterSpec>
struct GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > >
{
    typedef typename GetValue<TJournaledString>::Type Type;
};

template <typename TJournaledString, typename TJournaledStringIterSpec>
struct GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > const>
        : GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > > {};

template <typename TJournaledString, typename TJournaledStringIterSpec>
struct Reference<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > >
{
    typedef typename Reference<TJournaledString>::Type Type;
};

template <typename TJournaledString, typename TJournaledStringIterSpec>
struct Reference<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > const>
{
    typedef typename Reference<TJournaledString>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// For String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type
begin(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type TResult;
    TResult result;
    _initJournaledStringIterator(result, journalString);
    return result;
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type
begin(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type TResult;
    TResult result;
    _initJournaledStringIterator(result, journalString);
    return result;
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type
end(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journalString, Standard)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type TResult;
    TResult result;
    _initJournaledStringIteratorEnd(result, journalString);
    return result;
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type
end(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type TResult;
    TResult result;
    _initJournaledStringIteratorEnd(result, journalString);
    return result;
}

// For Iter<TJournaledString, JournaledStringIterSpec>

template <typename TJournaledString, typename TJournalSpec>
inline
void
_initJournaledStringIterator(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
                           TJournaledString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    iterator._journalEntriesIterator = begin(journalString._journalEntries);
    // Update iterators on the segment.
    _updateSegmentIterators(iterator);
}

template <typename TJournaledString, typename TJournalSpec>
inline
void
_initJournaledStringIteratorEnd(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
                              TJournaledString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    iterator._journalEntriesIterator = end(journalString._journalEntries);
}

// TODO(rmaerker): Rename to _updateSegmentIteratorsRight().
template <typename TJournaledString, typename TJournalSpec>
inline
void
_updateSegmentIterators(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (atEnd(iterator._journalEntriesIterator))
        return;
    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            iterator._hostSegmentBegin = begin(host(*iterator._journalStringPtr), Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._hostSegmentEnd = iterator._hostSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentHostIt = iterator._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            iterator._insertionBufferSegmentBegin = begin(iterator._journalStringPtr->_insertionBuffer, Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._insertionBufferSegmentEnd = iterator._insertionBufferSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentInsertionBufferIt = iterator._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
}

// ----------------------------------------------------------------------------
// Function _updateSegmentIteratorsLeft()
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TJournalSpec>
inline
void
_updateSegmentIteratorsLeft(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            iterator._hostSegmentBegin = begin(host(*iterator._journalStringPtr), Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._hostSegmentEnd = iterator._hostSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentHostIt = iterator._hostSegmentEnd - 1;
            break;
        case SOURCE_PATCH:
            iterator._insertionBufferSegmentBegin = begin(iterator._journalStringPtr->_insertionBuffer, Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._insertionBufferSegmentEnd = iterator._insertionBufferSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentInsertionBufferIt = iterator._insertionBufferSegmentEnd - 1;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline typename Reference<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type
value(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & me)
{
    SEQAN_CHECKPOINT;
    typedef typename Reference<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type TReference;
    TReference res(me);
    return res;
}

template <typename TJournaledString, typename TJournalSpec>
inline typename Reference<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const>::Type
value(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & me)
{
    SEQAN_CHECKPOINT;
    typedef typename Reference<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const>::Type TReference;
    TReference res(me);
    return res;
}

// assignValue

template <typename TJournaledString, typename TJournalSpec, typename TValue>
inline void
assignValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;
    assignValue(static_cast<TIterator const>(me), _value);
}

template <typename TJournaledString, typename TJournalSpec, typename TValue>
inline void
assignValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;
    typename Value<TIterator>::Type _temp_value = _value; //conversion
    assignValue(*iterator._journalStringPtr, position(iterator), _temp_value);
}

// moveValue

template <typename TJournaledString, typename TJournalSpec, typename TValue>
inline void
moveValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Copied from packed string. Actually, why no real move?
    assignValue(me, _value);
}

template <typename TJournaledString, typename TJournalSpec, typename TValue>
inline void
moveValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Copied from packed string. Actually, why no real move?
    assignValue(me, _value);
}

// valueConstruct

template <typename TJournaledString, typename TJournalSpec>
inline void
valueConstruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & /*it*/)
{
    // TODO(holtgrew): Intentionally left blank? Alphabet elements must be default-constructable.
}

template <typename TJournaledString, typename TJournalSpec, typename TParam>
inline void
valueConstruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & it,
               TParam const & param_)
{
    assignValue(it, param_);
}

template <typename TJournaledString, typename TJournalSpec, typename TParam>
inline void
valueConstruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & it,
               TParam & param_,
               Move const & /*tag*/)
{
    moveValue(it, param_);
}

// valueDestruct

template <typename TJournaledString, typename TJournalSpec>
inline void
valueDestruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & /*it*/)
{
    // TODO(holtgrew): Intentionally left blank? Copied from packed string, leads to problems with non-POD contents!
}

// position
// TODO(rmaerker): Implements an unexpected behavior. Should return the virtual position of the current iterator not the relative position within the current node.
template <typename TJournaledString, typename TJournalSpec>
inline typename Position<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const>::Type
position(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
{
    if (atEnd(iterator._journalEntriesIterator)) {
        return length(*iterator._journalStringPtr);
    }

    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            return iterator._currentHostIt - iterator._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            return iterator._currentInsertionBufferIt - iterator._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
            return 0;
    }
}

// ----------------------------------------------------------------------------
// Function _position()
// ----------------------------------------------------------------------------

// Returns the virtual position of the current iterator. Note, that the
// function position() should implement this behavior.
template <typename TJournaledString, typename TJournalSpec>
inline typename Position<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const>::Type
_position(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
{
    if (atEnd(iterator._journalEntriesIterator))
        return length(*iterator._journalStringPtr);

    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            return value(iterator._journalEntriesIterator).virtualPosition + iterator._currentHostIt - iterator._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            return value(iterator._journalEntriesIterator).virtualPosition + iterator._currentInsertionBufferIt - iterator._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
            return 0;
    }
}

// setPosition
template <typename TJournaledString, typename TJournalSpec, typename TPosition>
inline void
setPosition(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & /*me*/,
            TPosition /*pos_*/)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Implement me!
    SEQAN_ASSERT_FAIL("Set position...");
}

// getValue
template <typename TJournaledString, typename TJournalSpec>
inline
typename GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type
getValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return getValue(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return getValue(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type
getValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return getValue(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return getValue(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > &
operator++(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            ++iterator._currentHostIt;
            if (iterator._currentHostIt == iterator._hostSegmentEnd) {
                ++iterator._journalEntriesIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        case SOURCE_PATCH:
            ++iterator._currentInsertionBufferIt;
            if (iterator._currentInsertionBufferIt == iterator._insertionBufferSegmentEnd) {
                ++iterator._journalEntriesIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid Segment Source");
    }
    return iterator;
}

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
operator++(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator, int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > temp(iterator);
    ++iterator;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator--()                                               [prefix]
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > &
operator--(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    if (atEnd(iterator._journalEntriesIterator))
    {
        --iterator._journalEntriesIterator;
        _updateSegmentIteratorsLeft(iterator);
    }
    else
    {
        switch (value(iterator._journalEntriesIterator).segmentSource) {
            case SOURCE_ORIGINAL:
                if (iterator._currentHostIt == iterator._hostSegmentBegin) {
                    --iterator._journalEntriesIterator;
                    _updateSegmentIteratorsLeft(iterator);
                }
                else
                    --iterator._currentHostIt;
                break;
            case SOURCE_PATCH:
                if (iterator._currentInsertionBufferIt == iterator._insertionBufferSegmentBegin) {
                    --iterator._journalEntriesIterator;
                    _updateSegmentIteratorsLeft(iterator);
                }
                else
                    --iterator._currentInsertionBufferIt;
                break;
            default:
            {
                SEQAN_ASSERT_FAIL("Invalid segment source!");
            }
        }
    }
    return iterator;
}

// ----------------------------------------------------------------------------
// Function operator--()                                              [postfix]
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
operator--(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
            int /*postfix*/)
{
    Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > temp(iterator);
    --iterator;
    return temp;
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename Reference<TJournaledString>::Type
operator*(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename Reference<TJournaledString>::Type
operator*(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournaledString, typename TJournalSpec, typename TLen>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > &
operator+=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
           TLen len_)
{
    SEQAN_CHECKPOINT;

    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;

    // TODO(holtgrew): Handle case where len_ < 0?!
    SEQAN_ASSERT_GEQ(len_, static_cast<TLen>(0));
    size_t len = len_;

    // Handle bad case of len_ pointing at/behind end.
    if (position(iterator) + len_ >= length(*iterator._journalStringPtr)) {
        iterator = TIterator(end(*iterator._journalStringPtr));
        return iterator;
    }

    // Handle other case.
    typedef typename Size<TJournaledString>::Type TSize;
    while (len > 0) {
        TSize remaining;
        switch (value(iterator._journalEntriesIterator).segmentSource) {
            case SOURCE_ORIGINAL:
                remaining = iterator._hostSegmentEnd - iterator._currentHostIt;
                SEQAN_ASSERT_GT(remaining, 0u);
                if (len >= remaining) {
                    len -= remaining;
                    ++iterator._journalEntriesIterator;
                    _updateSegmentIterators(iterator);
                } else {
                    iterator._currentHostIt += len;
                    len = 0;
                }
                break;
            case SOURCE_PATCH:
                remaining = iterator._insertionBufferSegmentEnd - iterator._currentInsertionBufferIt;
                SEQAN_ASSERT_GT(remaining, 0u);
                if (len >= remaining) {
                    len -= remaining;
                    ++iterator._journalEntriesIterator;
                    _updateSegmentIterators(iterator);
                } else {
                    iterator._currentInsertionBufferIt += len;
                    len = 0;
                }
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid segment source!");
        }
    }
    return iterator;
}

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
operator+(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator,
          typename Size<TJournaledString>::Type const & len)
{
    SEQAN_CHECKPOINT;
    Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > temp(iterator);
    temp += len;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TJournalSpec, typename TLen>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > &
operator-=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
            TLen len_)
{
    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;
    typedef typename Position<TIterator>::Type TPosition;

    // TODO(holtgrew): Handle case where len_ < 0?!
    SEQAN_ASSERT_GEQ(len_, static_cast<TLen>(0));
    size_t len = len_;

    // Handle bad case of len_ pointing before begin.
    if (_position(iterator) <= static_cast<TPosition>(len_)) {
        iterator = TIterator(begin(*iterator._journalStringPtr));
        return iterator;
    }

    // Handle case when iterator is at positon at end
    if (atEnd(iterator._journalEntriesIterator))
    {
        --iterator._journalEntriesIterator;
        _updateSegmentIteratorsLeft(iterator);
        --len;
    }

    // Handle other case.
    typedef typename Size<TJournaledString>::Type TSize;
    while (len > 0) {
        TSize relNodePos;
        switch (value(iterator._journalEntriesIterator).segmentSource) {
            case SOURCE_ORIGINAL:
                relNodePos = iterator._currentHostIt - iterator._hostSegmentBegin;
                if (len > relNodePos) {
                    len -= (relNodePos + 1);
                    --iterator._journalEntriesIterator;
                    _updateSegmentIteratorsLeft(iterator);
                }
                else
                {
                    iterator._currentHostIt -= len;
                    len = 0;
                }
                break;
            case SOURCE_PATCH:
                relNodePos = iterator._currentInsertionBufferIt - iterator._insertionBufferSegmentBegin;
                if (len > relNodePos) {
                    len -= (relNodePos + 1);
                    --iterator._journalEntriesIterator;
                    _updateSegmentIteratorsLeft(iterator);
                } else {
                    iterator._currentInsertionBufferIt -= len;
                    len = 0;
                }
                break;
            default:
            {
                SEQAN_ASSERT_FAIL("Invalid segment source!");
            }
        }
    }
    return iterator;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
operator-(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator,
           typename Size<TJournaledString>::Type const & len)
{
    SEQAN_CHECKPOINT;
    Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > temp(iterator);
    temp -= len;
    return temp;
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename Difference<TJournaledString>::Type
operator-(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & it1,
          Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & it2)
{
    typedef typename Difference<TJournaledString>::Type TResult;

    // First, handle the cases where it1 or it2 are at the end.
    bool it1AtEnd = atEnd(it1._journalEntriesIterator);
    bool it2AtEnd = atEnd(it2._journalEntriesIterator);
    if (it1AtEnd && it2AtEnd) {
        return 0;
    } else if (it1AtEnd) {
        TResult len = length(*it1._journalStringPtr);
        TResult vPos = value(it2._journalEntriesIterator).virtualPosition;
        switch (value(it2._journalEntriesIterator).segmentSource) {
            case SOURCE_ORIGINAL:
                vPos += it2._currentHostIt - it2._hostSegmentBegin;
                break;
            case SOURCE_PATCH:
                vPos += it2._currentInsertionBufferIt - it2._insertionBufferSegmentBegin;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid segment source!");
        }
        SEQAN_ASSERT_LT(vPos, len);
        return len - vPos;
    } else if (it2AtEnd) {
        TResult len = length(*it1._journalStringPtr);
        TResult vPos = value(it1._journalEntriesIterator).virtualPosition;
        switch (value(it1._journalEntriesIterator).segmentSource) {
            case SOURCE_ORIGINAL:
                vPos += it1._currentHostIt - it1._hostSegmentBegin;
                break;
            case SOURCE_PATCH:
                vPos += it1._currentInsertionBufferIt - it1._insertionBufferSegmentBegin;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid segment source!");
        }
        SEQAN_ASSERT_LT(vPos, len);
        return vPos - len;
    }

    // Otherwise, we can simply subtract the virtual positions.
    TResult vPos1 = value(it1._journalEntriesIterator).virtualPosition;
    switch (value(it1._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            vPos1 += it1._currentHostIt - it1._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            vPos1 += it1._currentInsertionBufferIt - it1._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
    TResult vPos2 = value(it2._journalEntriesIterator).virtualPosition;
    switch (value(it2._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            vPos2 += it2._currentHostIt - it2._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            vPos2 += it2._currentInsertionBufferIt - it2._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
    return vPos1 - vPos2;
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator==(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    if (atEnd(a._journalEntriesIterator) && atEnd(b._journalEntriesIterator))
        return true;
    if (a._journalEntriesIterator != b._journalEntriesIterator)
        return false;
    if (value(a._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        if (a._currentHostIt != b._currentHostIt)
            return false;
    } else {
        SEQAN_ASSERT_EQ(value(a._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        if (a._currentInsertionBufferIt != b._currentInsertionBufferIt)
            return false;
    }
    return true;
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator==(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    typedef typename IterMakeConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type TConstIter;
    return static_cast<TConstIter>(a) == static_cast<TConstIter>(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator!=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    return !(a == b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator!=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return !(a == b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator<(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) < position(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator<(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) < position(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator<=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) <= position(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator<=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) <= position(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator>(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) > position(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator>(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) > position(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator>=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) >= position(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator>=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return position(a) >= position(b);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
