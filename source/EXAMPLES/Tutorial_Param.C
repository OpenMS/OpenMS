#include <OpenMS/DATASTRUCTURES/Param.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  Param param;
  
  param.setValue("file:name","test.xml");
  param.setValue("file:size(MB)",572.3);  
  param.setValue("file:data:min_int",0);
  param.setValue("file:data:max_int",16459);
  
  cout << "Name   : " << (String)(param.getValue("file:name")) << endl;
  cout << "Size   : " << (Real)(param.getValue("file:size(MB)")) << endl;
  cout << "Min int: " << (UInt)(param.getValue("file:data:min_int")) << endl;
  cout << "Max int: " << (UInt)(param.getValue("file:data:max_int")) << endl;
  
  return 0;
} //end of main
