/*
 * PackedMatrix.cpp
 *
 *  Created on: Feb 15, 2011
 *      Author: tomasoni
 */

#include "PackedMatrix.h"

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
  int row, col;
  for (row = 0; row < numRows(); row++) {
    for (col = 0; col < trhs.numRows(); col++) {
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
