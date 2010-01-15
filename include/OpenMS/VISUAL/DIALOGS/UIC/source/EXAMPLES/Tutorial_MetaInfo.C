#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  MetaInfoInterface info;
  
  //insert meta data
  info.setMetaValue("color",String("#ff0000"));
  info.setMetaValue("id",112131415);
  
  //access id by index
  UInt id_index = info.metaRegistry().getIndex("id");
  cout << "id   : " << (UInt)(info.getMetaValue(id_index)) << endl;
  //access color by name
  cout << "color: " << (String)(info.getMetaValue("color")) << endl;
  
  return 0;
} //end of main
