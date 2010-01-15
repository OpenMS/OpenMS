#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <iostream>

double someFunc(int, int*) { return 9.9; }

template < typename, int = 5 > struct WOW
{
  static float staticMemberFunc(float&) {return 0;}
  double memberFunc(const double &) {return 0;}
};

template < typename X, int I, unsigned J = 666, template < class, int > class WauWau = WOW > struct Oink {};

int main()
{
  using namespace OpenMS;
  using namespace std;

  cout << typeAsString(1) << endl;
  cout << typeAsString(2u) << endl;
  cout << typeAsString(3.) << endl;
  cout << typeAsString(4.f) << endl;
  cout << endl;
  cout << typeAsString<Int>() << endl;
  cout << typeAsString<PointerSizeUInt>() << endl;
  cout << endl;
  cout << typeAsString(Peak1D()) << endl;
  cout << typeAsString(DPeak<1>::Type()) << endl;
  cout << typeAsString(DPeak<1>::Type::PositionType()) << endl;
  cout << typeAsString<DPeak<1>::Type::CoordinateType>() << endl;
  cout << typeAsString<DPeak<1>::Type::IntensityType>() << endl;
  cout << endl;
  cout << typeAsString(&someFunc) << endl;
  cout << typeAsString<WOW<char const * const *** const & > >() << endl;
  cout << typeAsString<Oink<double,55> >() << endl;
  cout << typeAsString(&WOW<string,8>::staticMemberFunc) << endl;
  cout << typeAsString(&WOW<char,8>::memberFunc) << endl;

  return 0;
} // end of Tutorial_typeAsString.C
