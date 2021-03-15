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
// Journal entries implementation using an unbalanced binary search tree.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_H_
#define SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

// For iterator, needs to be declared here.
enum IterationDirection {
    DIRECTION_NULL,
    DIRECTION_DOWN_LEFT,
    DIRECTION_DOWN_RIGHT,
    DIRECTION_UP_LEFT,
    DIRECTION_UP_RIGHT
};


// Tag: Unbalanced tree.
struct UnbalancedTree {};


template <typename TNode, typename TTreeSpec>
class JournalEntries;


template <typename TCargo_>
class JournalEntries<TCargo_, UnbalancedTree>
{
  public:
    typedef TCargo_ TCargo;
    typedef JournalEntriesUnorderedTreeNode<TCargo> TNode;
    typedef JournalEntries<TCargo, UnbalancedTree> TJournalEntries;
    typedef typename Position<TNode>::Type TPos;
    typedef typename Size<TNode>::Type TSize;
    typedef UnbalancedTree TSpec;

    // Allocator for the nodes.
    Allocator<SinglePool<sizeof(TNode)> > _nodeAllocator;
    // Length of the underlying string.
    TSize _originalStringLength;
    // The root node.
    TNode * _root;

    JournalEntries()
            : _root(0)
    {
        SEQAN_CHECKPOINT;
    }

    JournalEntries(JournalEntries const & other)
            : _originalStringLength(other._originalStringLength)
    {
        _copyJournalEntriesNodes(_root, _nodeAllocator, other._root);
    }

    JournalEntries &
    operator=(JournalEntries const & other)

    {
        if (this != &other)
        {
            _originalStringLength = other._originalStringLength;
            _copyJournalEntriesNodes(_root, _nodeAllocator, other._root);
        }
        return *this;

    }

    JournalEntries &
    operator=(JournalEntries const & other) const

    {
        if (this != &other)
        {
            _originalStringLength = other._originalStringLength;
            _copyJournalEntriesNodes(_root, _nodeAllocator, other._root);
        }
        return *this;

    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TCargo>
struct Reference<JournalEntries<TCargo, UnbalancedTree> >
{
    typedef TCargo & Type;
};

template <typename TCargo>
struct Reference<JournalEntries<TCargo, UnbalancedTree> const>
{
    typedef TCargo const & Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TNode, typename TAllocator>
void
_copyJournalEntriesNodes(TNode * & target,
                      TAllocator & allocator,
                      TNode * const & source)
{
    SEQAN_CHECKPOINT;
    allocate(allocator, target, 1);
    target = new (target) TNode(*source);

    if (source->left != 0) {
        _copyJournalEntriesNodes(target->left, allocator, source->left);
        target->left->parent = target;
    }
    if (source->right != 0) {
        _copyJournalEntriesNodes(target->right, allocator, source->right);
        target->right->parent = target;
    }
}


template <typename TNode>
inline
bool checkVirtualPositionsRec(TNode * const & node, unsigned & virtualPosition)
{
    SEQAN_CHECKPOINT;
    if (node == 0)
        return true;
    bool res = true;
    if (node->left)
        res = res && checkVirtualPositionsRec(node->left, virtualPosition);
    if (cargo(*node).virtualPosition != virtualPosition)
        res = false;
    virtualPosition += cargo(*node).length;
    if (node->right)
        res = res && checkVirtualPositionsRec(node->right, virtualPosition);
    return res;
}

template <typename TNode>
inline
bool checkVirtualPositions(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    unsigned virtualPosition = 0;
    return checkVirtualPositionsRec(node, virtualPosition);
}

template <typename TNode>
inline
bool checkOrderRec(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    bool result = true;
    if (node->left != 0) {
        result = result && (cargo(*node).virtualPosition > cargo(*node->left).virtualPosition);
        result = result && checkOrderRec(node->left);
    }
    if (node->right != 0) {
        result = result && (cargo(*node).virtualPosition < cargo(*node->right).virtualPosition);
        result = result && checkOrderRec(node->right);
    }
    return result;
}

template <typename TNode>
inline
bool checkOrder(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    if (node == 0)
        return true;
    return checkOrderRec(node);
}

template <typename TNode>
inline
bool checkStructureRec(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    bool result = true;
    if (node->left != 0) {
        result = result && (node->left->parent == node);
        result = result && checkStructureRec(node->left);
    }
    if (node->right != 0) {
        result = result && (node->right->parent == node);
        result = result && checkStructureRec(node->right);
    }
    return result;
}

template <typename TNode>
inline
bool checkStructure(TNode * const & node)
{
    if (node == 0)
        return true;
    if (node->parent == 0)
        return checkStructureRec(node);
    if (node->parent->left != 0 && node->parent->right != 0) {
        if (!((node->parent->left == node) ^ (node->parent->right == node)))
            return false;
    } else if (node->parent->left == 0 && node->parent->right != 0) {
        if (node->parent->right != node)
            return false;
    } else if (node->parent->left != 0 && node->parent->right == 0) {
        if (node->parent->left != node)
            return false;
    } else {  // both == 0
        return false;
    }
    return checkStructureRec(node);
}

template <typename TStream, typename TNode>
inline
TStream &
operator<<(TStream & stream, JournalEntries<TNode, UnbalancedTree> const & tree)
{
    if (tree._root != 0)
        return stream << "JournalEntries(" << *tree._root << ")";
    else
        return stream << "JournalEntries()";
}

template <typename TCargo>
inline
void reinit(JournalEntries<TCargo, UnbalancedTree> & tree,
            typename Size<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type originalStringLength)
{
    SEQAN_CHECKPOINT;
    typedef typename JournalEntries<TCargo, UnbalancedTree>::TNode TNode;
    clear(tree._nodeAllocator);
    tree._originalStringLength = originalStringLength;
    TNode *tmp;
    allocate(tree._nodeAllocator, tmp, 1);
    tree._root = new (tmp) TNode(TCargo(SOURCE_ORIGINAL, 0, 0, 0, originalStringLength));
}


// Subtract delta from all nodes with virtual positions right of,
// respectively >= pos.  Note that this must not make the tree invalid.
template <typename TNode>
inline
void
_subtractFromVirtualPositionsRightOf(TNode * node,
                                     typename Position<TNode>::Type const & pos,
                                     typename Position<TNode>::Type const & delta)
{
    SEQAN_CHECKPOINT;

    if (cargo(*node).virtualPosition >= pos) {
        cargo(*node).virtualPosition -= delta;
        if (node->left != 0)
            _subtractFromVirtualPositionsRightOf(node->left, pos, delta);
        if (node->right != 0)
            _subtractFromVirtualPositionsRightOf(node->right, pos, delta);
    } else {  // node->virtualPosition < pos
        if (node->right != 0)
            _subtractFromVirtualPositionsRightOf(node->right, pos, delta);
    }
}


// Add delta to all nodes with virtual positions right of,
// respectively >= pos.  Note that this must not make the tree invalid.
template <typename TNode>
inline
void
_addToVirtualPositionsRightOf(TNode * node,
                              typename Position<TNode>::Type const & pos,
                              typename Position<TNode>::Type const & delta)
{
    SEQAN_CHECKPOINT;

    if (cargo(*node).virtualPosition >= pos) {
        cargo(*node).virtualPosition += delta;
        if (node->left != 0)
            _addToVirtualPositionsRightOf(node->left, pos, delta);
        if (node->right != 0)
            _addToVirtualPositionsRightOf(node->right, pos, delta);
    } else {  // node->virtualPosition < pos
        if (node->right != 0)
            _addToVirtualPositionsRightOf(node->right, pos, delta);
    }
}


template <typename TCargo>
inline
typename Iterator<JournalEntries<TCargo, UnbalancedTree> const, Standard>::Type
findInJournalEntries(JournalEntries<TCargo, UnbalancedTree> const & journalEntries,
                     typename Position<TCargo>::Type const & pos) {
    SEQAN_CHECKPOINT;
    typedef JournalEntries<TCargo, UnbalancedTree> TJournalEntries;
    typedef typename TJournalEntries::TNode TNode SEQAN_TYPEDEF_FOR_DEBUG;
    typedef typename Iterator<TJournalEntries, Standard>::Type TIterator;

    TIterator result;
    result._currentNode = journalEntries._root;
    result._iterationDirection = DIRECTION_DOWN_LEFT;

    while (true) {
        SEQAN_ASSERT_NEQ(result._currentNode, static_cast<TNode *>(0));

        if ((value(result).virtualPosition <= pos) && (value(result).virtualPosition + value(result).length > pos)) {
            break;
        } else if (value(result).virtualPosition + value(result).length <= pos) {
            if (!goRight(result))
                break;
        } else {  // pos < node->virtualPosition
            if (!goLeft(result))
                break;
        }
    }

    return result;
}


template <typename TCargo>
inline
TCargo const &
findJournalEntry(JournalEntries<TCargo, UnbalancedTree> const & journalEntries,
                 typename Position<TCargo>::Type const & pos) {
    SEQAN_CHECKPOINT;
    return *findInJournalEntries(journalEntries, pos);
}


template <typename TCargo>
inline
void recordErase(JournalEntries<TCargo, UnbalancedTree> & tree,
                 typename Position<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & pos,
                 typename Position<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & posEnd)
{
    SEQAN_CHECKPOINT;
    typedef JournalEntries<TCargo, UnbalancedTree> TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TIterator;
    typedef typename TJournalEntries::TNode TNode;
    typedef typename Position<TNode>::Type TPos;
    typedef typename Size<TNode>::Type TSize;

    SEQAN_ASSERT(checkStructure(tree._root));
    SEQAN_ASSERT(checkOrder(tree._root));
    SEQAN_ASSERT(checkVirtualPositions(tree._root));

//     std::cout << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    // Handle case of an empty journal tree.
    if (tree._root == 0) {
        SEQAN_ASSERT_EQ(pos, 0u);
        SEQAN_ASSERT_EQ(posEnd, 0u);
        return;
    }
    // Handle special case of removing all of the root node.
    if (tree._root->left == 0 && tree._root->right == 0 && pos == 0 && posEnd == cargo(*tree._root).length) {
        tree._root = 0;
        clear(tree._nodeAllocator);
        return;
    }

    // Find node with virtual position pos.
    TNode * node = 0;
    TNode * parent = 0;
    TIterator iter = findInJournalEntries(tree, pos);
    node = iter._currentNode;
    parent = node->parent;

    // The position to subtract values from right of.
    TPos subtractRightOf = pos;
    // Virtual begin and end position of the node.
    TPos nodeBegin = cargo(*node).virtualPosition;
    TPos nodeEnd = cargo(*node).virtualPosition + cargo(*node).length;

    // The simple cases are: A prefix or suffix but not all of the node are erased.
    // Otherwise, we have to delete or split the existing node.
    if (nodeBegin == pos && nodeEnd == posEnd) {
//         std::cout << "whole node" << std::endl;
        // The whole node is removed.  If there is <= one child, things are
        // simple, otherwise, we replace node with its right child and perform
        // a left-right traversal to find the leaf for which node's left child
        // can become the left child of.
        TNode * left = node->left;
        TNode * right = node->right;
        if (left == 0) {
            // Replace parent's pointer.
            if (parent == 0) {
                // Node was root.
                tree._root = node->right;
                tree._root->parent = 0;
            } else {
                // Node was not root.
                if (parent->left == node) {
                    parent->left = node->right;
                } else {
                    SEQAN_ASSERT_EQ(parent->right, node);
                    parent->right = node->right;
                }
                if (node->right != 0)
                    node->right->parent = parent;
            }
        } else {
            if (right == 0) { // left != 0 && right == 0
                if (parent == 0) {
                    // Node was root.
                    node->left->parent = 0;
                    tree._root = node->left;
                } else {
                    // Node was not root.
                    if (parent->left == node) {
                        node->left->parent = parent;
                        parent->left = node->left;
                    } else {
                        node->left->parent = parent;
                        parent->right = node->left;
                    }
                }
            } else {
                // left != 0 && right != 0
                if (parent == 0) {
                    // node is root
                    node->right->parent = 0;
                    tree._root = node->right;
                } else {
                    if (parent->left == node) {
                        node->right->parent = parent;
                        parent->left = node->right;
                    } else {
                        node->right->parent = parent;
                        parent->right = node->right;
                    }
                }
                TNode * tmp = node->right->left;
                node->right->left = node->left;
                node->left->parent = node->right;
                // Left-right traversal from node on.
                TNode * current = node->left;
                SEQAN_ASSERT_NEQ(current, static_cast<TNode *>(0));
                while (current->right != 0)
                    current = current->right;
                current->right = tmp;
                if (tmp != 0)
                    tmp->parent = current;
            }
        }
        // Actually deallocate node.
        deallocate(tree._nodeAllocator, node, 1);
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else if (nodeBegin == pos && nodeEnd > posEnd) {
//         std::cout << "true prefix" << std::endl;
        // A true prefix is removed, not the whole node.
        cargo(*node).length -= posEnd - pos;
        cargo(*node).physicalPosition += posEnd - pos;
        subtractRightOf = posEnd;  // No need to update this node!
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else if (nodeBegin < pos && nodeEnd == posEnd) {
//         std::cout << "true suffix" << std::endl;
        // A true suffix is removed, not the whole node.
        cargo(*node).length -= posEnd - pos;
        subtractRightOf = posEnd;  // No need to update this node!
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else if (nodeBegin < pos && nodeEnd > posEnd) {
//         std::cout << "true infix" << std::endl;
        // A true infix of the node is removed.  This node stores the prefix
        // and we do a right-left traversal to place the newly created node.
        TNode * tmp;
        allocate(tree._nodeAllocator, tmp, 1);
        // Perform index calculations for prefix and suffix node.
        TSize prefixLength = pos - cargo(*node).virtualPosition;
        TSize deletedInfixLength = posEnd - pos;
        TSize suffixLength = cargo(*node).length - prefixLength - deletedInfixLength;
        // Update prefix node.
        cargo(*node).length = prefixLength;
        // Construct suffix node.
        TNode * suffixNode = new (tmp) TNode(TCargo(cargo(*node).segmentSource,
                                                    cargo(*node).physicalPosition + prefixLength + deletedInfixLength,
                                                    cargo(*node).virtualPosition + prefixLength,
                                                    0,
                                                    suffixLength));
        // Insert node for suffix.
        if (node->right == 0) {
            node->right = suffixNode;
            suffixNode->parent = node;
        } else {
            TNode * current = node->right;
            while (current->left != 0) {
                current = current->left;
            }
            current->left = suffixNode;
            suffixNode->parent = current;
        }
        subtractRightOf = cargo(*node).virtualPosition + prefixLength + deletedInfixLength + 1;
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else {
//         std::cout << "spans more than one" << std::endl;
        // The deletion range spans more than this node.  We solve this by
        // calling this function recursively.  First, we delete a suffix of
        // this node (possibly all of it) and then we recursively remove the
        // remaining suffix of [pos, posEnd).
        TSize len = cargo(*node).length - (pos - cargo(*node).virtualPosition);
        recordErase(tree, pos, pos + len);
        recordErase(tree, pos, posEnd - len);
    }
//     unsigned nextId = mtRand();
//     journalTreeToDot(std::cerr, nextId, tree);
//     std::cout << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;
    SEQAN_ASSERT(checkStructure(tree._root));
    SEQAN_ASSERT(checkOrder(tree._root));
    SEQAN_ASSERT(checkVirtualPositions(tree._root));
}


template <typename TCargo>
inline
void recordInsertion(JournalEntries<TCargo, UnbalancedTree> & tree,
                     typename Position<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & virtualPos,
                     typename Position<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & physicalBeginPos,
                     typename Size<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & length)
{
    SEQAN_CHECKPOINT;
    typedef JournalEntries<TCargo, UnbalancedTree> TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TIterator;
    typedef typename TJournalEntries::TNode TNode;
    typedef typename Position<TNode>::Type TPos;

    SEQAN_ASSERT(checkStructure(tree._root));
    SEQAN_ASSERT(checkOrder(tree._root));
    SEQAN_ASSERT(checkVirtualPositions(tree._root));

    // Handle special case of empty tree.
    if (tree._root == 0) {
        SEQAN_ASSERT_EQ(virtualPos, 0u);
        if (length == 0)
            return;
        TNode * tmp;
        allocate(tree._nodeAllocator, tmp, 1);
        tree._root = new (tmp) TNode(TCargo(SOURCE_PATCH, physicalBeginPos, virtualPos, 0, length));
        return;
    }

    TNode * node;
    TNode * parent;
    TIterator iter = findInJournalEntries(tree, virtualPos);
    node = iter._currentNode;
    parent = node->parent;
    SEQAN_ASSERT_LEQ(cargo(*node).virtualPosition, virtualPos);

    if (cargo(*node).virtualPosition + cargo(*node).length > virtualPos) {
        // Found node that contains virtualPos.
        SEQAN_ASSERT_LEQ(cargo(*node).virtualPosition, virtualPos);
        if (cargo(*node).virtualPosition == virtualPos) {
            // Simple case: Insert left of current.
            TNode * tmp;
            allocate(tree._nodeAllocator, tmp, 1);
            TNode * insertNode = new (tmp) TNode(TCargo(SOURCE_PATCH, physicalBeginPos, virtualPos, 0, length));
            _addToVirtualPositionsRightOf(tree._root, virtualPos, length);
            insertNode->left = node->left;
            if (insertNode->left != 0)
                insertNode->left->parent = insertNode;
            node->left = insertNode;
            insertNode->parent = node;
        } else {
            // Harder case: Split current and insert new node.
            _addToVirtualPositionsRightOf(tree._root, cargo(*node).virtualPosition + cargo(*node).length, length);
            TPos offset = virtualPos - cargo(*node).virtualPosition;
            // Create right part of the node.
            TNode * tmp;
            allocate(tree._nodeAllocator, tmp, 1);
            TNode * splitNode = new (tmp) TNode(TCargo(cargo(*node).segmentSource, cargo(*node).physicalPosition + offset, cargo(*node).virtualPosition + offset + length, 0, cargo(*node).length - offset));
            // Current node becomes left part of current node.
            cargo(*node).length = offset;
            // Create insertion node.
            allocate(tree._nodeAllocator, tmp, 1);
            TNode * insertNode = new (tmp) TNode(TCargo(SOURCE_PATCH, physicalBeginPos, virtualPos, 0, length));
            // Insert into tree...
            insertNode->left = node;
            node->parent = insertNode;
            insertNode->right = splitNode;
            splitNode->parent = insertNode;
            splitNode->right = node->right;
            if (node->right != 0)
                node->right->parent = splitNode;
            node->right = 0;
            if (parent == 0) {
                // current is the root node.
                tree._root = insertNode;
                insertNode->parent = 0;
            } else {
                if (parent->left == node)
                    parent->left = insertNode;
                else
                    parent->right = insertNode;
                insertNode->parent = parent;
            }
        }
    } else {
        // Returned node with highest virtualPosition but we need to insert
        // right of it.
        SEQAN_ASSERT_EQ(node->right, static_cast<TNode *>(0));
        SEQAN_ASSERT_EQ(cargo(*node).virtualPosition + cargo(*node).length, virtualPos);
        TNode * tmp;
        allocate(tree._nodeAllocator, tmp, 1);
        TNode * insertNode = new (tmp) TNode(TCargo(SOURCE_PATCH, physicalBeginPos, virtualPos, 0, length));
        node->right = insertNode;
        insertNode->parent = node;
    }
//     unsigned nextId = mtRand();
//     journalTreeToDot(std::cerr, nextId, tree);
//     std::cout << tree << std::endl;
//     std::cout << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;
    SEQAN_ASSERT(checkStructure(tree._root));
    SEQAN_ASSERT(checkOrder(tree._root));
    SEQAN_ASSERT(checkVirtualPositions(tree._root));
}


// TODO(holtgrew): Remove?
/*
template <typename TNode>
inline
TNode const *
back(JournalEntries<TNode, UnbalancedTree> const & tree)
{
    SEQAN_XXXCHECKPOINT;
    TNode * current = tree._root;
    while (current->right != 0)
        current = current->right;
    return current;
}


template <typename TNode>
inline
TNode const *
front(JournalEntries<TNode, UnbalancedTree> const & tree)
{
    SEQAN_XXXCHECKPOINT;
    TNode * current = tree._root;
    while (current->left != 0)
        current = current->left;
    return current;
}
*/

//---------------------------------------------------
// function hostToVirtualPosition()
//---------------------------------------------------

template <typename TCargo, typename TPos>
inline
typename Position<TCargo>::Type
hostToVirtualPosition(JournalEntries<TCargo, UnbalancedTree> const & journalEntries,
                      TPos const & hostPos)
{
    typedef JournalEntries<TCargo, UnbalancedTree> TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TIterator;

    TIterator it = begin(journalEntries, Standard());
    SEQAN_ASSERT(it != end(journalEntries, Standard()));
    // TCargo tmp;
    while (it != end(journalEntries, Standard()))
    {
        // host position can only appear in orginal node
        if (value(it).segmentSource == SOURCE_ORIGINAL)
        {
            // found a potential hit - investigate!
            if (hostPos >= (TPos) value(it).physicalPosition)
            {
                // Case#1: Found searched position - report mapping to virtual position
                if (hostPos < (TPos) value(it).physicalPosition + (TPos) value(it).length)
                {
                    return value(it).virtualPosition + (hostPos - value(it).physicalPosition);
                }
                // Case#2a: Searched position appears either in next node (see Case#1) or
                // lies within a deletion (see Case#2b)
                // tmp = value(it);
            }
            else
            {   // Case#2b: Searched position lies within deletion - report virtual position of this node.
                return value(it).virtualPosition;
            }
        }
        ++it;
    }
    SEQAN_ASSERT_FAIL("Should never reach here!");
    return 0;
}



template <typename TStream, typename TNode>
void journalTreeToDotRec(TStream & stream, unsigned & nextId, TNode const & node)
{
    unsigned currentId = nextId;
    nextId += 1;
    stream << "  node_" << currentId << "[label=\"source=" << node.segmentSource << ", vpos=" << node.virtualPosition << ", phpos=" << node.physicalPosition << ", len=" << node.length << "\"]" << std::endl;
    if (node.left != 0) {
        stream << "  node_" << currentId << " -> node_" << nextId << "[label=\"L\"]" << std::endl;
        journalTreeToDotRec(stream, nextId, *(node.left));
        nextId += 1;
    }
    if (node.right != 0) {
        stream << "  node_" << currentId << " -> node_" << nextId << "[label=\"R\"]" << std::endl;
        journalTreeToDotRec(stream, nextId, *(node.right));
        nextId += 1;
    }
}


template <typename TStream, typename TNode>
void journalTreeToDot(TStream & stream, unsigned & nextId, JournalEntries<TNode, UnbalancedTree> const & journalTree)
{
    stream << "ROOTPTR" << nextId << " -> node_" << nextId << std::endl;
    journalTreeToDotRec(stream, nextId, *journalTree._root);
}


}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRIES_UNBALANCED_TREE_H_
