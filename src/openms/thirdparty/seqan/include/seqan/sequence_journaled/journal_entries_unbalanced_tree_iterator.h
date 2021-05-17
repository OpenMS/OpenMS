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
// Code for a SearchEntries<UnorderedTree> iterator.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_ITERATOR_H_
#define SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_ITERATOR_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

template <typename TJournalEntriesSpec>
struct JournalEntriesIterSpec;


template <typename TJournalEntries>
class Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> >
{
public:
    typedef Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > TIter;
    typedef typename TJournalEntries::TNode TNode;

    // The current node.
    TNode * _currentNode;
    // The direction we arrived from at this node.  The end node is
    // encoded as the root and _iterationDirection == DIRECTION_UP_RIGHT.
    IterationDirection _iterationDirection;

    Iter() : _currentNode(0), _iterationDirection(DIRECTION_NULL) {}

    Iter(TIter const & other)
            : _currentNode(other._currentNode),
              _iterationDirection(other._iterationDirection)
    { SEQAN_CHECKPOINT; }

    Iter(typename IterComplementConst<TIter>::Type const & other)
            : _currentNode(other._currentNode),
              _iterationDirection(other._iterationDirection)
    { SEQAN_CHECKPOINT; }

    explicit
    Iter(TJournalEntries & tree)
    {
        SEQAN_CHECKPOINT;
        _initJournalEntriesIterator(*this, tree);
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

template <typename TCargo>
struct Iterator<JournalEntries<TCargo, UnbalancedTree>, Standard>
{
    typedef Iter<JournalEntries<TCargo, UnbalancedTree>, JournalEntriesIterSpec<UnbalancedTree> > Type;
};

template <typename TCargo>
struct Iterator<JournalEntries<TCargo, UnbalancedTree> const, Standard>
{
    typedef Iter<JournalEntries<TCargo, UnbalancedTree> const, JournalEntriesIterSpec<UnbalancedTree> > Type;
};

template <typename TJournalEntries>
struct Value<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >
{
    typedef typename TJournalEntries::TCargo & Type;
};

template <typename TJournalEntries>
struct Value<Iter<TJournalEntries const, JournalEntriesIterSpec<UnbalancedTree> > >
{
    typedef typename TJournalEntries::TCargo const & Type;
};

template <typename TJournalEntries>
struct GetValue<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >
{
    typedef typename TJournalEntries::TCargo Type;
};

template <typename TJournalEntries>
struct GetValue<Iter<TJournalEntries const, JournalEntriesIterSpec<UnbalancedTree> > >
{
    typedef typename TJournalEntries::TCargo Type;
};

template <typename TJournalEntries>
struct Reference<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >
{
    typedef typename TJournalEntries::TCargo & Type;
};

template <typename TJournalEntries>
struct Reference<Iter<TJournalEntries const, JournalEntriesIterSpec<UnbalancedTree> > >
{
    typedef typename TJournalEntries::TCargo const & Type;
};

// ============================================================================
// Functions
// ============================================================================

// For JournalEntries<TNode, UnbalancedTree>

template <typename TNode>
inline
typename Iterator<JournalEntries<TNode, UnbalancedTree> const, Standard>::Type
begin(JournalEntries<TNode, UnbalancedTree> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return Iter<JournalEntries<TNode, UnbalancedTree> const, JournalEntriesIterSpec<UnbalancedTree> >(journalTree);
}

template <typename TNode>
inline
typename Iterator<JournalEntries<TNode, UnbalancedTree>, Standard>::Type
begin(JournalEntries<TNode, UnbalancedTree> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return Iter<JournalEntries<TNode, UnbalancedTree>, JournalEntriesIterSpec<UnbalancedTree> >(journalTree);
}

template <typename TNode>
inline
typename Iterator<JournalEntries<TNode, UnbalancedTree> const, Standard>::Type
end(JournalEntries<TNode, UnbalancedTree> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<JournalEntries<TNode, UnbalancedTree> const, Standard>::Type TIterator;
    TIterator result;
    result._currentNode = journalTree._root;
    result._iterationDirection = DIRECTION_UP_RIGHT;
    return result;
}

template <typename TNode>
inline
typename Iterator<JournalEntries<TNode, UnbalancedTree>, Standard>::Type
end(JournalEntries<TNode, UnbalancedTree> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef Iter<JournalEntries<TNode, UnbalancedTree>, JournalEntriesIterSpec<UnbalancedTree> > TIterator;
    TIterator result;
    result._currentNode = journalTree._root;
    result._iterationDirection = DIRECTION_UP_RIGHT;
    return result;
}

// For JournalEntriesIterator<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> >

template <typename TJournalEntries>
inline
void
_initJournalEntriesIterator(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator,
                         TJournalEntries & tree)
{
    SEQAN_CHECKPOINT;
    iterator._currentNode = tree._root;
    if (tree._root == 0) {
        iterator._iterationDirection = DIRECTION_UP_RIGHT;
    } else {
        iterator._iterationDirection = DIRECTION_DOWN_LEFT;
        while (goLeft(iterator))
            continue;  // Left-only traversal.
    }
}

// Initialize journal tree iterator to end, i.e. root & direction is up right.
template <typename TJournalEntries>
inline
void
_initJournalEntriesIteratorEnd(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator,
                            TJournalEntries & tree)
{
    SEQAN_CHECKPOINT;
    iterator._currentNode = tree._root;
    iterator._iterationDirection = DIRECTION_UP_RIGHT;
}

// TODO(holtgrew): Unused, remove?
/*
template <typename TJournalEntries>
inline
void
setValue(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator,
         typename Value<TJournalEntries>::Type const & value)
{
    SEQAN_XXXCHECKPOINT;
    setValue(*iterator._currentNode, value);
}
*/

template <typename TJournalEntries>
inline
typename Value<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >::Type
value(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_CHECKPOINT;
    return cargo(*iterator._currentNode);
}

template <typename TJournalEntries>
inline
typename Value<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >::Type
value(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return cargo(*iterator._currentNode);
}

// TODO(holtgrew): Unused, remove?
/*
template <typename TJournalEntries>
inline
typename Value<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >::Type
operator*(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_XXXCHECKPOINT;
    return value(iterator);
}
*/

// TODO(holtgrew): Unused, remove?
/*
template <typename TJournalEntries>
inline
typename Value<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >::Type
operator*(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & iterator)
{
    SEQAN_XXXCHECKPOINT;
    return value(iterator);
}
*/

template <typename TJournalEntries>
inline
bool
hasLeftChild(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_CHECKPOINT;
    return iterator._currentNode->left != 0;
}

template <typename TJournalEntries>
inline
bool
goLeft(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (!hasLeftChild(iterator))
        return false;
    iterator._iterationDirection = DIRECTION_DOWN_LEFT;
    iterator._currentNode = iterator._currentNode->left;
    return true;
}

template <typename TJournalEntries>
inline
bool
hasRightChild(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_CHECKPOINT;
    return iterator._currentNode->right != 0;
}

template <typename TJournalEntries>
inline
bool
goRight(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (!hasRightChild(iterator))
        return false;
    iterator._iterationDirection = DIRECTION_DOWN_RIGHT;
    iterator._currentNode = iterator._currentNode->right;
    return true;
}

template <typename TJournalEntries>
inline
bool
hasParent(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return iterator._currentNode->parent != 0;
}
    
template <typename TJournalEntries>
inline
bool
goUp(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    if (!hasParent(iterator)) {
        // Up from the root means we go to end, set direction to up-right.
        iterator._iterationDirection = DIRECTION_UP_RIGHT;
        return false;
    }
    if (iterator._currentNode->parent->left == iterator._currentNode)
        iterator._iterationDirection = DIRECTION_UP_LEFT;
    else
        iterator._iterationDirection = DIRECTION_UP_RIGHT;
    iterator._currentNode = iterator._currentNode->parent;
    return true;
}

template <typename TJournalEntries>
inline
bool
atEnd(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return (iterator._currentNode == 0) || (!hasParent(iterator) && (iterator._iterationDirection == DIRECTION_UP_RIGHT));
}

template <typename TJournalEntries>
inline
bool
atEnd(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    return (iterator._currentNode == 0) || (!hasParent(iterator) && (iterator._iterationDirection == DIRECTION_UP_RIGHT));
}

template <typename TJournalEntries>
inline
Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > &
operator++(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_CHECKPOINT;
    switch (iterator._iterationDirection) {
        case DIRECTION_DOWN_LEFT:
        case DIRECTION_DOWN_RIGHT:
            // Arrived here by going down (either right or left).  Next is
            // either right-left traversal (if node has a child to the right)
            // or going up if we do not have a child on the right.
            if (goRight(iterator)) {
                while (goLeft(iterator))
                    continue;
            } else {
                while (goUp(iterator) && iterator._iterationDirection == DIRECTION_UP_RIGHT)
                    continue;
            }
            break;
        case DIRECTION_UP_LEFT:
            // We came up from our left child.  Next is either right-left
            // traversal (if node has a child to the right) or going up if
            // we do not have a child on the left.
            if (goRight(iterator)) {
                while (goLeft(iterator))
                    continue;
            } else {
                while (goUp(iterator) && iterator._iterationDirection == DIRECTION_UP_RIGHT)
                    continue;
            }
            break;
        case DIRECTION_UP_RIGHT:
            // We came here from our right child.  Next is up until we got
            // to the current node from it's left child.
            while (goUp(iterator) && iterator._iterationDirection == DIRECTION_UP_RIGHT)
                continue;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid iteration direction.");
    }
        
    return iterator;    
}

template <typename TJournalEntries>
inline
Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> >
operator++(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator,
           int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > temp(iterator);
    ++iterator;
    return temp;
}

template <typename TJournalEntries>
inline
Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > &
operator--(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator)
{
    SEQAN_CHECKPOINT;
    typedef typename TJournalEntries::TNode TNode;
    
    // If we are at the end, go to the last vertex.
    if (atEnd(iterator)) {
        while (goRight(iterator))
            continue;
        return iterator;
    }
    // Otherwise, simply go one step left.
    
    // The following comments are from the point of perspective of having
    // performed an increment in the previous step.
    TNode const * parent;
    switch (iterator._iterationDirection) {
        case DIRECTION_DOWN_RIGHT:
        case DIRECTION_DOWN_LEFT:
            // Arrived here in a right-left traversal.  Reverse this by
            // going up until we went up from a right child.
            while (goUp(iterator) && iterator._iterationDirection != DIRECTION_UP_RIGHT)
                continue;
            // Next, update iteration direction such that going up and then
            // down the same edge would lead to current node.
            parent = iterator._currentNode->parent;
            if (parent == 0) {
                iterator._iterationDirection = DIRECTION_DOWN_LEFT;
            } else {
                if (parent->left == iterator._currentNode) {
                    iterator._iterationDirection = DIRECTION_DOWN_LEFT;
                } else {
                    SEQAN_ASSERT_EQ(parent->right, iterator._currentNode);
                    iterator._iterationDirection = DIRECTION_DOWN_RIGHT;
                }
            }
            break;
        case DIRECTION_UP_LEFT:
            // We came up from our left child.  Next in going to predecessor
            // is going back to the left child.
            goLeft(iterator);
            break;
        case DIRECTION_UP_RIGHT:
            SEQAN_ASSERT_FAIL("Should probably not happen.");
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid iteration direction.");
    }
        
    return iterator;    
}

template <typename TJournalEntries>
inline
Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> >
operator--(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > & iterator,
           int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > temp(iterator);
    --iterator;
    return temp;
}

template <typename TJournalEntries>
inline
bool
operator==(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & a,
           Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode == b._currentNode) && (a._iterationDirection == b._iterationDirection);
}

template <typename TJournalEntries>
inline
bool
operator==(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & a,
           typename IterComplementConst<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode == b._currentNode) && (a._iterationDirection == b._iterationDirection);
}

template <typename TJournalEntries>
inline
bool
operator!=(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & a,
           Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode != b._currentNode) || (a._iterationDirection != b._iterationDirection);
}

template <typename TJournalEntries>
inline
bool
operator!=(Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > const & a,
           typename IterComplementConst<Iter<TJournalEntries, JournalEntriesIterSpec<UnbalancedTree> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode != b._currentNode) || (a._iterationDirection != b._iterationDirection);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_ITERATOR_H_
