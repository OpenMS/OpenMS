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
// Code for the node of JournalEntries<UnorderedTree>, i.e. a node of an
// unbalanced binary search tree.  This could easily be extended to a Red Black
// Tree with a flag (which could actually be stored in the lowest bits of one
// of the pointers).  However, since std::set already implements balanced
// trees, it should be sufficient to use a STL implementation of std::set for
// balanced binary search trees.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_NODE_H_
#define SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_NODE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

template <typename TCargo>
struct JournalEntriesUnorderedTreeNode
{
    // Left child.
    JournalEntriesUnorderedTreeNode * left;
    // Right child.
    JournalEntriesUnorderedTreeNode * right;
    // Parent, 0 for root.
    JournalEntriesUnorderedTreeNode * parent;
    // The actual payload:  The journal entry.
    TCargo cargo;
    
    JournalEntriesUnorderedTreeNode() : left(0), right(0), parent(0) {}

    JournalEntriesUnorderedTreeNode(TCargo const & _cargo)
            : left(0),
              right(0),
              parent(0),
              cargo(_cargo) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// TODO(holtgrew): Rename value to cargo?

template <typename T>
struct Cargo;

template <typename TCargo>
struct Cargo<JournalEntriesUnorderedTreeNode<TCargo> >
{
    typedef TCargo Type;
};

template <typename TCargo>
struct Cargo<JournalEntriesUnorderedTreeNode<TCargo> const>
{
    typedef TCargo const Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TCargo>
TStream & operator<<(TStream & stream, JournalEntriesUnorderedTreeNode<TCargo> const & node)
{
    stream << "JournalEntriesUnorderedTreeNode(add=" << &node
           << ", cargo=" << node.cargo
           << ", parent=" << node.parent
           << ", left=";
    if (node.left)
        stream << *node.left;
    else
        stream << "NULL";
    stream << ", right=";
    if (node.right)
        stream << *node.right;
    else
        stream << "NULL";
    return stream << ")";
}

template <typename TCargo>
typename Cargo<JournalEntriesUnorderedTreeNode<TCargo> >::Type &
cargo(JournalEntriesUnorderedTreeNode<TCargo> & node)
{
    SEQAN_CHECKPOINT;
    return node.cargo;
}

// TODO(holtgrew): Unused, remove?
/*
template <typename TCargo>
typename Cargo<JournalEntriesUnorderedTreeNode<TCargo> const>::Type &
cargo(JournalEntriesUnorderedTreeNode<TCargo> const & node)
{
    SEQAN_XXXCHECKPOINT;
    return node.cargo;
}

template <typename TCargo>
typename Cargo<JournalEntriesUnorderedTreeNode<TCargo> const>::Type
getCargo(JournalEntriesUnorderedTreeNode<TCargo> const & node)
{
    SEQAN_XXXCHECKPOINT;
    return node.cargo;
}
*/

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_NODE_H_
