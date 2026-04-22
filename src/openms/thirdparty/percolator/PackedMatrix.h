/*
 * PackedMatrix.h
 *
 *  Created on: Feb 15, 2011
 *      Author: tomasoni
 */

#ifndef _PackedMatrix_H
#define _PackedMatrix_H

#include "PackedVector.h"
#include "Array.h"

class PackedMatrix
{
  public:
    PackedMatrix() {
      rows = Array<PackedVector>();
      cols = 0;
    }
    /**build a PackedMatrix without initializing the elements in the underlying
     * PackedVectors
     */
    PackedMatrix(int r, int c) {
      rows = Array<PackedVector>(r, PackedVector(0));
      cols = c;
    }
    virtual ~PackedMatrix(){};

    PackedMatrix packedTranspose(const PackedMatrix& mat);
    PackedMatrix packedAdd(const PackedMatrix& rhs);
    PackedMatrix packedMultiply(double scale);
    PackedMatrix packedMultiply(const PackedMatrix & rhs);
    PackedVector packedMultiply(const PackedVector & rhs);
    static PackedMatrix packedDiagonalMatrix(const Vector & v);
    void displayMatrix() const;

    PackedVector & operator [](int k){
      return rows[k];
    }
    const PackedVector & operator [](int k) const {
      return rows[k];
    }
    PackedVector & operator [](std::size_t k){
      return rows[k];
    }
    const PackedVector & operator [](std::size_t k) const {
      return rows[k];
    }
    int numRows() const {
      return static_cast<int>(rows.size());
    }
    int numCols() const {
      return cols;
    }

  protected:
    int cols;
    Array<PackedVector> rows;
};

#endif
