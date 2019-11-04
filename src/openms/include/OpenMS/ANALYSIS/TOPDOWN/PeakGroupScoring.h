//
// Created by Kyowon Jeong on 10/8/19.
//

#ifndef OPENMS_HOST_PEAKGROUPSCORING_H
#define OPENMS_HOST_PEAKGROUPSCORING_H

#include <Eigen/Dense>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
namespace OpenMS
{
  class OPENMS_DLLAPI PeakGroupScoring{
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalcularedAveragine PrecalcularedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    PeakGroupScoring(std::vector<PeakGroup> &peakGroups, Parameter &param);

    /// default destructor
    ~PeakGroupScoring();


    static double getChargeFitScore(double *perChargeIntensity, int range);

    static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                           double *perIsotopeIntensities,
                                                           int perIsotopeIntensitiesSize,
                                                           int &offset,
                                                           FLASHDeconvHelperStructs::PrecalcularedAveragine &avg);

    std::vector<PeakGroup> & scoreAndFilterPeakGroups(FLASHDeconvHelperStructs::PrecalcularedAveragine &avg);

  protected:
    std::vector<PeakGroup> &peakGroups;
    Parameter &param;

    void removeOverlappingPeakGroups();

    std::vector<int> updatePerChargeIsotopeIntensity(
        double **intensityGrid,
        double **intensityGrid2,
        double *perIsotopeIntensity,
        double *perChargeIntensity,
        PeakGroup &pg);

    static bool checkChargeDistribution(double *perChargeIntensity, int range, int threshold);


    static double getCosine(std::vector<double> &a, std::vector<double> &b, int off = 0);

    static double getCosine(double *a, double *b, int size);

    static double getCorrelation(double *a,
                                 int &aStart,
                                 int &aEnd,
                                 IsotopeDistribution &b,
                                 int &bSize,
                                 int offset);

    static double getCosine(double *a,
                            int &aStart,
                            int &aEnd,
                            IsotopeDistribution &b,
                            int &bSize,
                            double &bNorm,
                            int offset);
  };




}
#endif //OPENMS_HOST_PEAKGROUPSCORING_H
