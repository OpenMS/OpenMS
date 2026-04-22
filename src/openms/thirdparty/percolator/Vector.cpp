// Written by Oliver Serang 2009
// see license for more information

#include "Vector.h"

Numerical Vector::sparseChecker(1e-15);
Numerical Vector::comparator(1e-8);

#define PACK
//#define REPACK

const Vector & Vector::operator =(const Array<double> & simpleVector)
{
  values = simpleVector;
  createNonzeroIndices();
  return *this;
}

Vector::Vector(const Array<double> & simpleVector, bool pack)
{
  values = simpleVector;
  createNonzeroIndices();
}

const Vector & Vector::operator +=(const Vector & rhs)
{
  Set zeroedIndices;
  for (Set::Iterator iter = rhs.beginNonzero(); iter != rhs.endNonzero(); iter++)
    {
      values[ *iter ] += rhs[ *iter ];
      if ( sparseChecker.isZero( values[ *iter ] ) )
	zeroedIndices.add( *iter );
    }

  #ifdef REPACK
  nonzeroIndices = (nonzeroIndices | rhs.nonzeroIndices).without(zeroedIndices);
  #else
  nonzeroIndices = (nonzeroIndices | rhs.nonzeroIndices);
  #endif

  return *this;
}

const Vector & Vector::operator -=(const Vector & rhs)
{
  Set zeroedIndices;
  for (Set::Iterator iter = rhs.beginNonzero(); iter != rhs.endNonzero(); iter++)
    {
      values[ *iter ] -= rhs[ *iter ];
    }

  #ifdef REPACK
  nonzeroIndices = (nonzeroIndices | rhs.nonzeroIndices).without(zeroedIndices);
  #else
  nonzeroIndices = (nonzeroIndices | rhs.nonzeroIndices);
  #endif

  return *this;
}

const Vector & Vector::addEqScaled(double coef, const Vector & rhs)
{
  Set zeroedIndices;
  for (Set::Iterator iter = rhs.beginNonzero(); iter != rhs.endNonzero(); iter++)
    {
      values[ *iter ] += coef * rhs[ *iter ];
      if ( sparseChecker.isZero( values[ *iter ] ) )
	zeroedIndices.add( *iter );
    }

  #ifdef REPACK
  nonzeroIndices = (nonzeroIndices | rhs.nonzeroIndices).without(zeroedIndices);
  #else
  nonzeroIndices = (nonzeroIndices | rhs.nonzeroIndices);
  #endif

  return *this;
}

Vector Vector::operator [](const Set & rhs) const
{
  Array<double> subValues = values[ rhs ];
  
  Vector v = subValues;
  return v;
}

Vector operator +(const Vector & lhs, const Vector & rhs)
{
  Vector result = lhs;
  result += rhs;
  return result;
}

Vector operator -(const Vector & lhs, const Vector & rhs)
{
  Vector result = lhs;
  result -= rhs;
  return result;
}

Vector Vector::operator -() const
{
  Vector result = *this;

  for (Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++)
    {
      result.values[ *iter ] = -result[ *iter ];
    }

  return result;
}

void Vector::createNonzeroIndices()
{
#ifdef PACK
  nonzeroIndices = Set();
  
  for (int k=0; k<size(); k++)
    {
      if ( sparseChecker.isNonzero( values[k] ) )
	{
	  nonzeroIndices.add( k );
	}
    }
#else
  nonzeroIndices = Set::FullSet(0 , size()-1 );
#endif
}

void Vector::trimNonzeroIndices()
{
#ifdef REPACK
  Set newNonzeroIndices;
  for (Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++)
    {
      if ( sparseChecker.isNonzero( values[ *iter ] ) )
	{
	  newNonzeroIndices.add( *iter );
	}
    }
  nonzeroIndices = newNonzeroIndices;
#endif
}

Array<double> Vector::unpack() const
{
  return values;
}

Vector operator *(double val, const Vector & rhs)
{
  Vector result = rhs;
  result *= val;
  return result;
}

Vector operator /(const Vector & rhs, double val)
{
  Vector result = rhs;
  result /= val;
  return result;
}

Set Vector::operator <(double val) const
{
  Set result;
  
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      if ( comparator.isNeg( values[ *iter ] - val ) )
	result.add( *iter );
    }

  if ( comparator.isNeg( 0.0 - val ) )
    {
      // add the zero value indices
      result |= Set::FullSet(0, size()-1 ).without( nonzeroIndices );
    }

  return result;
}

Set Vector::operator >(double val) const
{
  Set result;
  
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      if ( comparator.isPos( values[ *iter ] - val ) )
	result.add( *iter );
    }

  if ( comparator.isPos( 0.0 - val ) )
    {
      // add the zero value indices
      result |= Set::FullSet(0, size()-1 ).without( nonzeroIndices );
    }

  return result;
}

Set Vector::operator <=(double val) const
{
  Set result;
  
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      if ( comparator.isNonpos( values[ *iter ] - val ) )
	result.add( *iter );
    }

  if ( comparator.isNonpos( 0.0 - val ) )
    {
      // add the zero value indices
      result |= Set::FullSet(0, size()-1 ).without( nonzeroIndices );
    }

  return result;
}

Set Vector::operator >=(double val) const
{
  Set result;
  
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      if ( comparator.isNonneg( values[ *iter ] - val ) )
	result.add( *iter );
    }

  if ( comparator.isNonneg( 0.0 - val ) )
    {
      // add the zero value indices
      result |= Set::FullSet(0, size()-1 ).without( nonzeroIndices );
    }

  return result;
}


bool Vector::prec(double val) const
{
  return ( (*this) >= val ).isEmpty();
}

bool Vector::succ(double val) const
{
  return ( (*this) <= val ).isEmpty();
}

bool Vector::precEq(double val) const
{
  return ( (*this) > val ).isEmpty();
}

bool Vector::succEq(double val) const
{
  return ( (*this) < val ).isEmpty();
}

Set Vector::operator ==(double val) const
{
  Set result;
  
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      if ( comparator.isZero( values[ *iter ] - val ) )
	result.add( *iter );
    }

  if ( comparator.isZero( 0.0 - val ) )
    {
      // add the zero value indices
      result |= Set::FullSet(0, size()-1 ).without( nonzeroIndices );
    }

  return result;
}

Set Vector::operator !=(double val) const
{
  Set result;
  
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      if ( comparator.isNonzero( values[ *iter ] - val ) )
	result.add( *iter );
    }

  if ( comparator.isNonzero( 0.0 - val ) )
    {
      // add the zero value indices
      result |= Set::FullSet(0, size()-1 ).without( nonzeroIndices );
    }

  return result;
}

double operator *(const Vector & lhs, const Vector & rhs)
{
  if ( lhs.size() != rhs.size() )
    throw Vector::DimensionException();
          
  if ( lhs.numberEntries() < rhs.numberEntries() )
    return rhs * lhs;

  // assume that rhs is the sparser of the two vectors

  double tot = 0;

  for ( Set::Iterator iterRhs = rhs.beginNonzero(); iterRhs != rhs.endNonzero(); iterRhs++ )
    {
      tot += lhs[ *iterRhs ] * rhs[ *iterRhs ];
    }

  return tot;
}

Vector operator /(const Vector & lhs, const Vector & rhs)
{
  if ( lhs.size() != rhs.size() )
    throw Vector::DimensionException();

  Vector result( lhs.size() );

  for ( Set::Iterator iter = lhs.beginNonzero(); iter != lhs.endNonzero(); iter++ )
    {
      result.addElement( *iter, lhs[ *iter ] / rhs[ *iter ]);
    }

  return result;
}

ostream & operator <<(ostream & os, const Vector & rhs)
{
  os << "( " << rhs.size() << " : " << rhs.nonzeroIndices << " ; " << rhs.values[ rhs.nonzeroIndices ] << " )";
  return os;
}

// Commented away this statement 240312
/*
istream & operator >>(istream & is, Vector & rhs)
{
  char c;
  
  is >> c;
  if ( c != '(' )
    throw Vector::FormatException();

  int n;
  is >> n;
  rhs = Vector(n);

  is >> c;
  if ( c != ':' )
    throw Vector::FormatException();

  is >> rhs.nonzeroIndices;

  is >> c;
  if ( c != ';' )
    throw Vector::FormatException();

  Array<double> values;
  is >> values;

  is >> c;
  if ( c != ')' )
    throw Vector::FormatException();

  int counter = 0;
  for (Set::Iterator iter = rhs.beginNonzero(); iter != rhs.endNonzero(); iter++, counter++)
    {
      rhs.values[ *iter ] = values[ counter ];
    }

  return is;
}
*/

double operator *(const Array<double> & lhs, const Array<double> & rhs)
{
  double tot = 0;
  for (int k=0; k<lhs.size(); k++)
    {
      tot += lhs[k] * rhs[k];
    }

  return tot;
}

const Array<double> & operator *=(Array<double> & lhs, double val)
{
  for (int k=0; k<lhs.size(); k++)
    lhs[k] *= val;
  return lhs;
}

double norm(const Array<double> & vec)
{
  return sqrt(vec * vec);
}

void Vector::add(double val)
{
  resize( size()+1 );
  addElement( size() - 1, val);
}

void Vector::addElement(int ind, double val)
{
  if ( sparseChecker.isNonzero(val) )
    {
      nonzeroIndices.add( ind );
      values[ ind ] = val;
    }
}

void Vector::resize(int newSize)
{
  values.resize( newSize );
}

const Vector & Vector::operator *=(double val)
{
  for (Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++)
    {
      values[ *iter ] *= val;
    }

  return *this;
}

const Vector & Vector::operator /=(double val)
{
  return (*this) *= (1/val);
}

Vector Vector::normalized() const
{
  return (*this) / norm();
}

double Vector::min() const
{
  double currentMin = Numerical::inf();

  Set zeroIndices = Set::FullSet(0, size()-1 ).without( nonzeroIndices );
 
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      currentMin = ::min( currentMin, values[ *iter ] );     
    }
 
  if ( ! zeroIndices.isEmpty() )
    {
      currentMin = ::min( currentMin, 0.0 );
    }

  return currentMin;
}

double Vector::max() const
{
  double currentMax = -Numerical::inf();

  Set zeroIndices = Set::FullSet(0, size()-1 ).without( nonzeroIndices );
 
  for ( Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++ )
    {
      currentMax = ::max( currentMax, values[ *iter ] );     
    }
 
  if ( ! zeroIndices.isEmpty() )
    {
      currentMax = ::max( currentMax, 0.0 );
    }

  return currentMax;
}

Set Vector::argmin() const
{
  return (*this) == min();
}

Set Vector::argmax() const
{
  return (*this) == max();
}

// performance note: all of these could be better
double Vector::sum() const
{
  return ( Vector( Array<double>( size(), 1.0 ) ) * (*this) );
}

double Vector::prod() const
{
  double res = 1.0;

  for (Set::Iterator iter = beginNonzero(); iter != endNonzero(); iter++)
    {
      res *= (*this)[ *iter ];
    }

  return res;
}

double Vector::average() const
{
  return sum() / size();
}

double Vector::variance() const
{
  double mean = average();
  
  Array<double> deviations( size() );
  for (int k=0; k<size(); k++)
    {
      deviations[k] = pow( (*this)[k] - mean , 2.0 );
    }
  
  return Vector(deviations).sum() / ( size() - 1 );
}

Array<double> operator -(const Array<double> & rhs)
{
  Array<double> result = rhs;

  for (int k=0; k<result.size(); k++)
    {
      result[k] = -result[k];
    }
  
  return result;
}

void Vector::displayVector() const
{
  Array<double> u = unpack();

  for (int k=0; k<u.size(); k++)
    {
      cout << u[k];
      if ( k != u.size() - 1 )
	{
	  cout << "\t";
	}
    }
  cout << endl;
}

double diffNormSquared(const Vector & u, const Vector & v)
{
  double res = 0.0;

  Set s = u.nonzeroIndices & v.nonzeroIndices;
  
  for (int i=0; i<s.size(); i++)
    {
      int ind = s[i];
      res += pow(u[ind] - v[ind], 2.0);
    }
  return sqrt(res);
}
