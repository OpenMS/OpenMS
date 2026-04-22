// Written by Oliver Serang 2009
// see license for more information

#ifndef _Vector_H
#define _Vector_H
#include <cmath>
#include "Set.h"
#include "Numerical.h"
#include "Array.h"

class Vector
{
public:
  Vector()
    {
    }
 explicit Vector(int N) :
  values(N, 0.0)
    {
    }
  explicit Vector(int N, double val) :
  values(N, val)
    {
      createNonzeroIndices();
    }
  Vector(const Array<double> & rhs)
    {
      *this = rhs;
    }
  Vector(const std::vector<double> &rhs)
   {
      *this = Array<double>(rhs);
   }
    
  Vector(const Array<double> & rhs, bool repack);

  Array<double> unpack() const;

  const Vector & operator =(const Array<double> & rhs);

  const Vector & operator +=(const Vector & rhs);
  const Vector & operator -=(const Vector & rhs);
  const Vector & addEqScaled(double coef, const Vector & rhs);
  const Vector & operator *=(double val);
  const Vector & operator /=(double val);

  const double & operator [](int k) const
  {
    return values[k];
  }

  const double & operator [](std::size_t k) const
  {
    return values[k];
  }

  
  Set::Iterator beginNonzero() const
    {
      // decide whether to use full or not...
      return nonzeroIndices.begin();
    }

  Set::Iterator endNonzero() const
    {
      return nonzeroIndices.end();
    }

  // vector scaling
  friend Vector operator *(double val, const Vector & rhs);
  friend Vector operator /(const Vector & rhs, double val);

  void displayVector() const;
  
  Vector operator -() const;

  int size() const
  {
    return static_cast<int>(values.size());
  }

  int numberEntries() const
  {
    return static_cast<int>(nonzeroIndices.size());
  }
 
  friend Vector operator /(const Vector & lhs, const Vector & rhs);

  // inner product
  friend double operator *(const Vector & lhs, const Vector & rhs);


  // Set-wise indexing
  Vector operator [](const Set & rhs) const;

  friend ostream & operator <<(ostream & os, const Vector & rhs);
  //friend istream & operator >>(istream & is, Vector & rhs);
  
  double norm() const
  {
    return sqrt( (*this) * (*this) );
  }

  Set operator <(double val) const;
  Set operator >(double val) const;
  Set operator <=(double val) const;
  Set operator >=(double val) const;
  Set operator ==(double val) const;
  Set operator !=(double val) const;

  bool prec(double val) const;
  bool succ(double val) const;
  bool precEq(double val) const;
  bool succEq(double val) const;

  double max() const;
  double min() const;

  Set argmax() const;
  Set argmin() const;

  void verifyVectorIntegrity() const;

  void resize(int n);

  void add(double val);
  void addElement(int ind, double val);

  Vector normalized() const;
  double sum() const;
  double prod() const;
  double average() const;
  double variance() const;

  class DivisionByZeroException {};
  class FormatException {};
  class VectorIntegrityException {};
  class DimensionException {};
  class VectorInverseException {};
  class SmallEmbedException {};
  class EmptyVectorException {};
  class VectorIndexOfOfBoundsException {};

  static Numerical comparator;
  static Numerical sparseChecker;

  void trimNonzeroIndices();
  void createNonzeroIndices();
  
  friend double diffNormSquared(const Vector & u, const Vector & v);

//protected:

  Set nonzeroIndices;
  Array<double> values;
  
};

Vector operator +(const Vector & lhs, const Vector & rhs);
Vector operator -(const Vector & lhs, const Vector & rhs);

double operator *(const Array<double> & lhs, const Array<double> & rhs);
const Array<double> & operator *=(Array<double> & lhs, double val);
double norm(const Array<double> & vec);

Array<double> operator -(const Array<double> & rhs);


#endif

