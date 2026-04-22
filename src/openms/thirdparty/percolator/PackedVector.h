/**
 * @file PackedVector.h
 * @author Mattia Tomasoni
 * @section DESCRIPTION
 * The class PackedVector augments the class Vector with algebraic operations
 * on packed Vectors.
 */

#ifndef _PackedVector_H
#define _PackedVector_H

#define PACK
#define ALLOW_ZERO
#include "Vector.h"

class PackedVector : public Vector
{
  public:
    /**
     * Constructs an empty PackedVector
     */
    PackedVector() : Vector() {}
    /**
     * Constructs a PackedVector with elements of value 0 and creates the
     * corresponding indices
     *
     * @param N number of elements
     */
    explicit PackedVector(int N) : Vector(N) {
      createIndices();
    }
    /**
     * Constructs a PackedVector and initializes the corresponding indices
     *
     * @param N number of elements
     * @param val value assigned to the elements
     */
    explicit PackedVector(int N, double val) : Vector(N, val) {}

    /**
     * Replaces an existing element
     *
     * @param ind index of the element to be replaced
     * @param val new value to be assigned to the element
     */
    void packedReplace(int ind, double val);

    /**
     * Adds an element and updates the indexes
     *
     * @ param ind index of the new element (to be added)
     * @ param val value to be assigned to the new element
     */
    void packedAddElement(int ind, double val);

    /**
     * Returns the ith index among nonzeroIndices
     *
     * @param i position in the nonzeroIndices to be returned
     * @return ith index among nonzeroIndices
     */
    int index(int i) const;

    /**
     * Returns the position of index of value x in nonzeroIndices
     *
     * @param x value of the sought index
     * @return position in the nonzeroIndices. -1 if x is no present
     */
    int find(int x) const;

    /**
     * Swaps two elements in the values Array
     *
     * @param ind1 position of the first element in values
     * @param ind2 position of the second element in values
     */
    void swapElements(int ind1, int ind2);

    /**
     * Returns a non-packed (sparse) representation of the object where missing
     * indices in nonzeroIndices have been added together with a 0-value
     * element in the values Array
     */
    PackedVector makeSparse() const;

    PackedVector packedProd(double val);
    PackedVector packedDiv(double val);
    PackedVector packedSubtract(const PackedVector & rhs) const;
    PackedVector packedAdd(const PackedVector & rhs) const;
    double packedDotProd(const PackedVector& rhs) const;

    void displayVector() const;

    friend double packedNorm(const PackedVector& vec);

  //private:
    void createIndices();
};

double packedNorm(const PackedVector& vec);

/**
 * Returns true when two PackedVectors have the same values associated with the
 * same indices
 */
bool operator == (const PackedVector & lhs, const PackedVector & rhs);

#endif
