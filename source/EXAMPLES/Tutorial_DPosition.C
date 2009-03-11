#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <iostream>

using namespace OpenMS;

Int main()
{
  DPosition<2> pos;
  pos[0] = 8.15;
  pos[1] = 47.11;

  for (Size i = 0; i < DPosition<2>::DIMENSION; ++i)
  {
    std::cout << "Dimension " << i << ": "<< pos[i] << std::endl;
  }  

  return 0;
} //end of main
