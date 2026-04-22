/*
 * PackedVector.cpp
 *
 *  Created on: Feb 15, 2011
 *      Author: tomasoni
 */

#include "Numerical.h"
#include "PackedVector.h"
#include "assert.h"

PackedVector PackedVector::packedProd(double val) {
  for (int k = 0; k<size(); k++)
    packedReplace(k, values[k]*val);
  return *this;
}

PackedVector PackedVector::packedDiv(double val) {
  assert(sparseChecker.isNonzero(val));
  for (int k = 0; k < size(); k++)
    packedReplace(k, values[k]/val);
  return *this;
}

PackedVector PackedVector::packedSubtract(const PackedVector & rhs) const{
  PackedVector result(0);
  int kL, kR;
  for (kL = 0, kR = 0; kL < numberEntries() && kR < rhs.numberEntries();) {
    int iL = index(kL);
    int iR = rhs.index(kR);
    if (iL < iR) {
      result.packedAddElement(iL, values[kL]);
      kL++;
    } else if (iL > iR) {
      result.packedAddElement(iR, -rhs[kR]);
      kR++;
    }else {
      double res = values[kL] - rhs[kR];
      if ( sparseChecker.isNonzero(res) ) {
        result.packedAddElement(iL, res);
      }
      kL++;
      kR++;
    }
  }
  // add anything left over
  if (kL < numberEntries()) {
    for (; kL < numberEntries(); kL++) {
      result.packedAddElement(index(kL), values[kL]);
    }
  }
  if (kR < rhs.numberEntries()) {
    for (; kR < rhs.numberEntries(); kR++) {
      result.packedAddElement(rhs.index(kR), -rhs[kR]);
    }
  }
  assert(result.size()==result.numberEntries());
  return result;
}

PackedVector PackedVector::packedAdd(const PackedVector & rhs) const{
  PackedVector neg = rhs;
  neg.packedProd(-1);
  return packedSubtract(neg);
}

double PackedVector::packedDotProd(const PackedVector& rhs) const{
  double tot = 0;
  int kL, kR;
  for (kL = 0, kR = 0; kL < numberEntries() && kR < rhs.numberEntries();) {
    int iL = index(kL);
    int iR = rhs.index(kR);
    if (iL < iR) {
      kL++;
    } else if (iL > iR) {
      kR++;
    } else {
      // equality case
      tot += values[find(iL)] * rhs[rhs.find(iR)];
      kL++;
      kR++;
    }
  }
  return tot;
}

int PackedVector::index(int i) const {
  return nonzeroIndices[i];
}

int PackedVector::find(int x) const {
  return nonzeroIndices.find(x);
}

void PackedVector::packedReplace(int ind, double val)
{
  assert(ind<size());
  assert(size()==numberEntries());
  values[ ind ] = val;
}

void PackedVector::packedAddElement(int ind, double val)
{
  assert(nonzeroIndices.find(ind) == -1);
  nonzeroIndices.add(ind);
  values.add(val);
  assert(size()==numberEntries());
}

void PackedVector::swapElements(int ind1, int ind2){
  double oldInd1 = values[ind1];
  values[ind1] = values[ind2];
  values[ind2] = oldInd1;
}

void PackedVector::createIndices(){
  nonzeroIndices = Set();
  for (int k=0; k<size(); k++){
    nonzeroIndices.add( k );
  }
  assert(size()==numberEntries());
}

PackedVector PackedVector::makeSparse() const{
  PackedVector sparse;
  int count(0);
  int ind(0);
  for(;ind<numberEntries();){
    int sparseInd = index(ind);
    if(sparseInd == count){
      sparse.packedAddElement(count, values[ind]);
      ind++; count++;
    } else {
      sparse.packedAddElement(count, 0);
      count++;
    }
  }
  assert(size()==numberEntries());
  return sparse;
}

void PackedVector::displayVector() const {
  assert(size()==numberEntries());
  for (int k=0; k<size(); k++) {
    cerr << nonzeroIndices[k] << ":" << values[k];
    if (k!=size()-1) cerr << "\t";
  }
  cerr << endl;
}

double packedNorm(const PackedVector & vec){
  return sqrt(vec.packedDotProd(vec));
}

bool operator == (const PackedVector & lhs, const PackedVector & rhs){
  assert(lhs.size()==lhs.numberEntries());
  assert(rhs.size()==rhs.numberEntries());
  bool different = false;
  if(lhs.numberEntries()!=rhs.numberEntries()) different = true;
  Set::Iterator it = lhs.beginNonzero();
  while (!different && it!=lhs.endNonzero()){
    double lhs_value = lhs[lhs.find(*it)];
    double rhs_value = rhs[rhs.find(*it)];
    if(lhs_value == rhs_value) it++;
    else different = true;
  }
  return !different;
}
