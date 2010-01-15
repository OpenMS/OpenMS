#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  const ElementDB* db = ElementDB::getInstance();

  Element carbon = *db->getElement("Carbon"); // .getResidue("C") would also be ok

  cout   << carbon.getName() << " " 
        << carbon.getSymbol() << " " 
        << carbon.getMonoWeight() << " "
        << carbon.getAverageWeight() << endl;

  return 0;
} //end of main
