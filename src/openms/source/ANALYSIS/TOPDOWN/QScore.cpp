
#include "OpenMS/ANALYSIS/TOPDOWN/QScore.h"
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

namespace OpenMS
{

  double QScore::getQScore(PeakGroup *pg, double intensity, int charge){
    double score;
   if (pg == nullptr)
    { // all zero
      score = -1.3923 * log10(intensity + 1) + 8.3952;
    }
    else
    {
      score =
          -0.367 * log10(pg->perChargeInfo[charge][0] + 1e-3)
          - 1.4339 * log10(intensity + 1)
          + 1.3297 * log10(pg->perChargeInfo[charge][2] + 1)
          -2.7199 * pg->perChargeInfo[charge][1]
          - 0.8385 * log10(pg->totalSNR + 1e-3)
          + 4.29 * pg->isotopeCosineScore
          - 0.4913 * pg->chargeCosineScore
          - 0.4561 * log10(pg->intensity + 1)
          + 1.424;
      /*
auto score =
        -0.367 * log10(pg.perChargeSNR[charge] + 1e-3)
        - 1.4339 * log10(preChargeMaxIntensity[j] + 1)
        + 1.3297 * log10(pg.perChargeSumInt[charge] + 1)
        -2.7199 * pg.perChargeICos[charge]
        - 0.8385 * log10(pg.totalSNR + 1e-3)
        + 4.29 * pg.isotopeCosineScore
        - 0.4913 * pg.chargeCosineScore
        - 0.4561 * log10(pg.intensity + 1)
        + 1.424;

          -0.5495 * log10(pg->perChargeSNR[peak->charge] + 1e-3)
              - 0.3263 * log10(peak->intensity + 1)
              - 0.9337 * log10(pg->totalSNR + 1e-3)
              + 2.8442 * pg->isotopeCosineScore
              - 0.9721 * pg->chargeCosineScore
              - 0.1885 * log10(pg->intensity + 1)
              + 1.9703;*/
    }

    return -score;
  }

  //double QScore::getQScore(PeakGroup *pg, Precursor peak){
  //  return getQScore(pg, peak.getIntensity(), peak.getCharge());
  //}


}
