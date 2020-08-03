//
// Created by Kyowon Jeong on 10/8/19.
//

#pragma once

#include <Eigen/Dense>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

namespace OpenMS
{
  class OPENMS_DLLAPI PeakGroupScoring
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    //typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
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
                                                           FLASHDeconvHelperStructs::PrecalculatedAveragine &avg);

    std::vector<PeakGroup> &scoreAndFilterPeakGroups(unsigned int &msLevel,
                                                     FLASHDeconvHelperStructs::PrecalculatedAveragine &avg);

  protected:
    std::vector<PeakGroup> &peakGroups;
    Parameter &param;


    void removeOverlappingPeakGroups(double tol);
    void removeHarmonicPeakGroups(double tol);
    std::vector<int> updatePerChargeIsotopeIntensity(
        //        double **intensityGrid,
        //        double **intensityGrid2,
        double *perIsotopeIntensity,
        double *perChargeIntensity,
        PeakGroup &pg);

    void filterPeakGroupsByIsotopeCosine(int currentMaxMassCount);
    void filterPeakGroupsByQScore(int currentMaxMassCount);

    void filterPeakGroupsByIntensity(int currentMaxMassCount);

    double getPredictionScore(PeakGroup &pg, int charge); //

    //static double getAvgMassPpmError(PeakGroup &pg);


    static bool checkChargeDistribution(double *perChargeIntensity, int range, int threshold);


    static double getCosine(std::vector<double> &a, std::vector<double> &b, int off = 0);

    static double getCosine(const double *a, double *b, Size size);

    static double getCorrelation(const double *a,
                                 int &aStart,
                                 int &aEnd,
                                 IsotopeDistribution &b,
                                 int &bSize,
                                 int offset);

    static double getCosine(const double *a,
                            int &aStart,
                            int &aEnd,
                            IsotopeDistribution &b,
                            int &bSize,
                            double &bNorm,
                            int offset);
  public:

  };


}
