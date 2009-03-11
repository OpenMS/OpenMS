#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <iostream>

using namespace OpenMS;

Int main()
{
  DRange<2> range;  
  range.setMin( DPosition<2>(2.0, 3.0) );
  range.setMax( DPosition<2>(1.0, 5.0) );

  for (UInt i = 0; i < DRange<2>::DIMENSION; ++i)
  {
    std::cout << "min " << i << ": "<< range.min()[i] << std::endl;
    std::cout << "max " << i << ": "<< range.max()[i] << std::endl;
  }  

  return 0;
} //end of main
