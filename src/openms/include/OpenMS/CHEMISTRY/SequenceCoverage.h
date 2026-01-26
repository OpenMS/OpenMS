#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>

namespace OpenMS
{


  class OPENMS_DLLAPI SequenceCoverage
  {
  public:
 
    static double getCoverage(
      const AASequence& protein,
      const std::vector<AASequence>& peptides
    );
  };

}
