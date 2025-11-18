from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/BinaryTreeNode.h>" namespace "OpenMS":

    cdef cppclass BinaryTreeNode:
        BinaryTreeNode(Size left_child, Size right_child, float distance) except + nogil
        BinaryTreeNode(BinaryTreeNode&) except + nogil
        Size left_child
        Size right_child
        float distance
