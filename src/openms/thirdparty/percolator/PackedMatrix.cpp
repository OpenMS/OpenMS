/*
 * PackedMatrix.cpp
 *
 *  Created on: Feb 15, 2011
 *      Author: tomasoni
 */

#include "PackedMatrix.h"

#include <algorithm>
#include <vector>
namespace OpenMS { namespace Internal { namespace Percolator {

PackedMatrix PackedMatrix::packedTranspose(const PackedMatrix& mat) {
  PackedMatrix res = PackedMatrix(mat.numCols(), mat.numRows());
  int i, j;
  for (i = 0; i < mat.numRows(); ++i) {
    for (j = 0; j < mat[i].numberEntries(); ++j) {
      res[mat[i].index(j)].packedAddElement(i, mat[i][j]);
    }
  }
  return res;
}

PackedMatrix PackedMatrix::packedAdd(const PackedMatrix& rhs) {
  PackedMatrix res = PackedMatrix(numRows(),numCols());
  int row;
  for (row = 0; row < numRows(); row++) {
    res[row] = rows[row].packedAdd(rhs[row]);
  }
  return res;
}

PackedMatrix PackedMatrix::packedMultiply(double scale)
{
  PackedMatrix result = *this;
  for (int k=0; k<result.numRows(); k++)
  {
    result[k].packedProd(scale);
  }
  return result;
}

PackedMatrix PackedMatrix::packedMultiply(const PackedMatrix & rhs){
  PackedMatrix res = PackedMatrix(numRows(),rhs.numCols());
  PackedMatrix trhs = trhs.packedTranspose(rhs);
  // OpenMS: per row, only the columns that share an index with it (the nonzero entries of the
  // rows of rhs at the row's indices), in increasing order: the dot product with any other
  // column is 0, which is not stored. The products of the (banded) matrices of BaseSpline take
  // linear instead of quadratic time; the result is the same.
  std::vector<int> cols;
  int row;
  for (row = 0; row < numRows(); row++) {
    cols.clear();
    for (int k = 0; k < rows[row].numberEntries(); ++k) {
      const PackedVector& rhsRow = rhs[rows[row].index(k)];
      for (int m = 0; m < rhsRow.numberEntries(); ++m) {
        cols.push_back(rhsRow.index(m));
      }
    }
    std::sort(cols.begin(), cols.end());
    cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
    for (int col : cols) {
      double prod = rows[row].packedDotProd(trhs[col]);
      if (Vector::sparseChecker.isNonzero(prod)) {
        res[row].packedAddElement(col, prod);
      }
    }
  }
  return res;
}

PackedVector PackedMatrix::packedMultiply(const PackedVector & rhs){
  PackedVector res;
  int row;
  for (row = 0; row < numRows(); row++) {
    double prod = rows[row].packedDotProd(rhs);
    res.packedAddElement(row, prod);
  }
  return res;
}

PackedMatrix PackedMatrix::packedDiagonalMatrix(const Vector & v ){
  PackedMatrix result( v.size() , v.size());
  for (int k=0; k<v.size(); k++)
    result.rows[k].packedAddElement( k, v[k] );
  return result;
}

/** displaying indexes as well as values for packed matrixes
 */
void PackedMatrix::displayMatrix() const{
  cerr << "{";
  for (int row=0; row<numRows(); row++){
    cerr << "{ ";
    for (int k = 0; k < rows[row].numberEntries(); k++) {
      cerr << rows[row].index(k) << ":" << rows[row][k];
      if (k != rows[row].numberEntries() - 1) {
        cerr << ", ";
      }
    }
    cerr << " }";

    if (row != numRows() - 1) {
      cerr << "," << endl;
    }
  }
  cerr << "}"<<endl;
}

}}}  // namespace OpenMS::Internal::Percolator
