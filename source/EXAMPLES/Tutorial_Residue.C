#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  const ResidueDB* res_db = ResidueDB::getInstance();

  Residue lys = *res_db->getResidue("Lysine"); // .getResidue("K") would also be ok

  cout   << lys.getName() << " " 
        << lys.getThreeLetterCode() << " " 
        << lys.getOneLetterCode() << " "
        << lys.getAverageWeight() << endl;

  return 0;
} //end of main
